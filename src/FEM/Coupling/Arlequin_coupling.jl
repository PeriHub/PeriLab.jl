# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Arlequin_coupling
include("../Element_formulation/Lagrange_element.jl")
include("../FEM_routines.jl")
using .Lagrange_element:
    define_lagrangian_grid_space, get_recursive_lagrange_shape_functions
using .FEM_routines: get_polynomial_degree
include("../../Support/Helpers.jl")
using .Helpers:
    find_active_nodes,
    find_point_in_polygon,
    find_point_in_hexagon,
    get_ring,
    get_hexagon,
    create_centroid_and_search_radius,
    get_nearest_neighbors
using LinearAlgebra


export coupling_name
export init_coupling_model
export compute_coupling
function coupling_name()
    return "Arlequin"
end

function init_coupling_model(datamanager::Module, nodes, fe_params::Dict)
    @info "Coupling $(coupling_name()) is active."
    dof = datamanager.get_dof()
    p = get_polynomial_degree(fe_params, dof)
    if !haskey(fe_params["Coupling"], "PD Weight")
        fe_params["Coupling"]["PD Weight"] = 0.5
        @info "Weight for coupling is set to 0.5."
    end

    el_topology = datamanager.get_field("FE el_topology")
    coordinates = datamanager.get_field("Coordinates")

    pd_nodes = datamanager.get_field("PD Nodes")
    fe_nodes = datamanager.get_field("FE Nodes")
    pd_nodes = find_active_nodes(fe_nodes, pd_nodes, nodes, false)
    if haskey(fe_params["Coupling"], "Coupling Block")
        @info "Pre defined coupling region is block $(fe_params["Coupling"]["Coupling Block"])"
        pd_nodes = get_coupling_zone(
            pd_nodes,
            datamanager.get_field("Block_Id"),
            fe_params["Coupling"]["Coupling Block"],
        )
    end
    lumped_mass = datamanager.get_field("Lumped Mass Matrix")

    # TODO memory efficiency; create a mapping field and allocate only the memory of length coupling nodes
    coupling_matrix = datamanager.create_constant_node_field(
        "Coupling Matrix",
        Float64,
        "Matrix",
        (prod(p .+ 1) + 1),
    )

    if any(x -> x > 1, p)
        @error "Coupling is supported only for linear elements yet."
        return nothing
    end
    if !haskey(fe_params["Coupling"], "Kappa")
        fe_params["Coupling"]["Kappa"] = 1.0
        @info "Coupling stiffness is set to 1."
    end
    kappa = fe_params["Coupling"]["Kappa"]
    weight_coefficient = fe_params["Coupling"]["PD Weight"]
    @info "Coupling stiffness kappa: $kappa"
    # Assign Element ids
    @info "Find Coupling Pairs"

    coupling_dict = find_point_in_elements(coordinates, el_topology, pd_nodes, dof)
    # point should be inside element -> for Arlequin; more stable
    rho = datamanager.get_field("Density")
    @info "Create Coupling Matrix"
    for (coupling_node, coupling_element) in pairs(coupling_dict)
        coupling_matrix[coupling_node, :, :] = compute_coupling_matrix(
            coordinates,
            el_topology,
            coupling_node,
            coupling_element,
            kappa,
            p,
            dof,
        )
        rho[coupling_node] *= (1 - weight_coefficient)
    end
    datamanager.set_coupling_dict(coupling_dict)
    element_list = collect(values(coupling_dict))
    unique_fe_coupling_nodes = unique(el_topology[element_list, :])
    lumped_mass[unique_fe_coupling_nodes] .*= weight_coefficient

    datamanager.set_coupling_fe_nodes(unique_fe_coupling_nodes)
    return datamanager
end


function get_coupling_zone(pd_nodes, blocks, ref_block)
    return [pd_nodes[i] for i in eachindex(blocks) if blocks[i] == ref_block]
end
function topo_closed_loop(pn)
    mapping = Int64[]

    for ip = 1:pn[1]+1
        push!(mapping, ip)
    end
    for ip = 2:pn[2]+1
        push!(mapping, ip * (pn[1] + 1))
    end

    for ip = 2:pn[1]+1
        push!(mapping, (pn[1] + 1) * (pn[2] + 1) - ip + 1)
    end

    for ip = 1:pn[2]-1
        push!(mapping, (pn[1] + 1) * pn[2] + 1 - (pn[1] + 1) * ip)
    end

    return mapping
end
##
# TODO
# only one point PD -> one Element; Check
# weights are not multiplied for FE -> lumped and coupling; check for internal -> kuerzt sich raus?
# 0.25 shape function check # TODO 0.25  -> not in FEM
##

function compute_coupling(datamanager::Module, fem_params::Dict)
    dof = datamanager.get_dof()
    el_topology = datamanager.get_field("FE el_topology")
    force_densities = datamanager.get_field("Force Densities", "NP1")
    displacements = datamanager.get_field("Displacements", "NP1")
    volume = datamanager.get_field("Volume")
    coupling_matrix = datamanager.get_field("Coupling Matrix")
    # EQ(13) in 10.1002/pamm.202400021 is fulfilled , also for the mass part.
    # F_i/rho = a ?! TODO check this assumption
    weight_coefficient = fem_params["Coupling"]["PD Weight"]
    coupling_dict = datamanager.get_coupling_dict()
    unique_fe_coupling_nodes = datamanager.get_coupling_fe_nodes()
    force_densities[unique_fe_coupling_nodes, :] .*= weight_coefficient

    # TODO weights are not correct here see EQ(1)
    for (coupling_node, coupling_element) in pairs(coupling_dict)

        topo = vcat(coupling_node, el_topology[coupling_element, :])
        force_densities[coupling_node, :] .*= (1 - weight_coefficient)
        vol = vcat(volume[coupling_node], volume[el_topology[coupling_element, :]])

        force_densities[topo, :] -=
            coupling_matrix[coupling_node, :, :] * displacements[topo, :] ./ vol

        #force_densities[topo, :] +=

        #force_densities[coupling_node, :] += force_temp[1, :] ./ volume[coupling_node]
        #force_densities[el_topology[coupling_element, :], :] += force_temp[2:5, :] ./ volume[coupling_node]^2
        ## TODO det(jacobi)
    end

    return datamanager

end


function compute_coupling_matrix(
    coordinates,
    el_topology,
    point_val,
    element_number,
    kappa,
    p,
    dof,
)
    # only one point per call
    # Point coordinates
    x_point = coordinates[point_val, 1]
    y_point = coordinates[point_val, 2]
    #
    ### Right now only for linear elements
    x1 = coordinates[el_topology[element_number, 1], 1]
    y1 = coordinates[el_topology[element_number, 1], 2]
    x2 = coordinates[el_topology[element_number, 2], 1]
    y2 = coordinates[el_topology[element_number, 2], 2]
    x3 = coordinates[el_topology[element_number, 4], 1]
    y3 = coordinates[el_topology[element_number, 4], 2]
    x4 = coordinates[el_topology[element_number, 3], 1]
    y4 = coordinates[el_topology[element_number, 3], 2]
    #### this is not correct, because diag2 is dx
    dx = sqrt((x2 - x1)^2 + (y2 - y1)^2)

    # Compute diagonal lengths
    diag1 = sqrt((x3 - x1)^2 + (y3 - y1)^2)
    diag2 = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    diag3 = sqrt((x4 - x1)^2 + (y4 - y1)^2)

    if diag1 > dx
        ksi = 2 * (x_point - (x3 + x1) / 2) / dx
        eta = 2 * (y_point - (y3 + y1) / 2) / dx
    elseif diag2 > dx
        ksi = 2 * (x_point - (x2 + x1) / 2) / dx
        eta = 2 * (y_point - (y2 + y1) / 2) / dx
    elseif diag3 > dx
        ksi = 2 * (x_point - (x4 + x1) / 2) / dx
        eta = 2 * (y_point - (y4 + y1) / 2) / dx
    else
        # Fallback to center or default values
        ksi = 0.0
        eta = 0.0
    end

    xi = define_lagrangian_grid_space(dof, p)
    Nxi = get_recursive_lagrange_shape_functions(xi[1, :], ksi, p[1])
    Neta = get_recursive_lagrange_shape_functions(xi[2, :], eta, p[2])
    # shape functions of PD point in local coord
    # notation from Lagrange_element.jl
    N1p = Nxi[1] * Neta[1] #(1 - ksi) * (1 - eta)
    N2p = Nxi[2] * Neta[1] #(1 + ksi) * (1 - eta)
    N3p = Nxi[1] * Neta[2] #(1 - ksi) * (1 + eta)
    N4p = Nxi[2] * Neta[2] #(1 + ksi) * (1 + eta)

    Np = [N1p N2p N3p N4p]
    local_coupl_matrix = kappa * [1 -Np; -Np' Np'*Np]

    return local_coupl_matrix
end


function find_point_in_elements(coordinates, el_topology, points_to_check, dof)
    if dof == 2
        fu = get_ring
        inside_check = find_point_in_polygon
    else
        fu = get_hexagon
        inside_check = find_point_in_hexagon
    end
    el_centroid, search_radius =
        create_centroid_and_search_radius(coordinates, el_topology, dof, fu)
    near_points = fill(Vector{Int64}([]), length(search_radius))
    near_points = get_nearest_neighbors(
        1:length(search_radius),
        dof,
        el_centroid,
        points_to_check,
        search_radius,
        near_points,
    )
    # nearest neighbors -> centroid ID is equal to the element

    return find_point_in_element(el_topology, near_points, coordinates, fu, Dict())
end


## TODO does only work if the PD point is not on a FE - point, line or surface (3D)



end
