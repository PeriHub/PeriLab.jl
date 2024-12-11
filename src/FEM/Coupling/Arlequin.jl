# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Arlequin_coupling
include("../Element_formulation/lagrange_element.jl")
include("../FEM_routines.jl")
using .Lagrange_element:
    define_lagrangian_grid_space, get_recursive_lagrange_shape_functions
using .FEM_routines: get_polynomial_degree
using .Helpers: find_active_nodes
using LinearAlgebra
function coupling_name()
    return "Arlequin"
end

function init_coupling_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    complete_params::Dict,
)
    dof = datamanager.get_dof()
    p = get_polynomial_degree(complete_params, dof)
    if !haskey(complete_params, "Weight")
        complete_params["PD Weight"] = 0.5
        @info "Weight for coupling is set to 0.5."
    end
    topology = datamanager.get_field("FE Topology")
    coordinates = datamanager.get_field("Coordinates")

    pd_nodes = datamanager.get_field("PD Nodes")
    fe_nodes = datamanager.get_field("FE Nodes")
    pd_nodes = find_active_nodes(fe_nodes, pd_nodes, nodes)
    datamanager.create_constant_node_field("Number of Element Neighbors", Int64, 1, 1)
    # Find number of element neighbors
    # TODO size correct?

    coupling_matrix = datamanager.create_constant_node_field(
        "Coupling Matrix",
        Float64,
        "Matrix",
        prod(p .+ 1) + 1,
        prod(p .+ 1) + 1,
    )
    # TODO to specify over params?
    kappa = 1
    # Assign Element ids
    coupling_dict = find_point_in_elements(coordinates, topology, pd_nodes)
    for (coupling_node, coupling_element) in pairs(coupling_dict)
        coupling_matrix[coupling_node, :, :] = compute_coupling_matrix(
            coordinates,
            topology,
            coupling_node,
            coupling_element,
            kappa,
            p,
            dof,
        )
    end
    datamanager.set_coupling_dict(coupling_dict)
    return datamanager
end

function compute_coupling(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    fem_params::Dict,
)

    topology = datamanager.get_field("FE Topology")
    force_densities = datamanager.get_field("Force Densities", "NP1")
    displacements = datamanager.get_field("Displacements", "NP1")
    weight_coefficient = fem_params["PD Weight"]
    coupling_matrix = datamanager.get_field("Coupling Matrix")
    coupling_dict = datamanager.get_coupling_dict()

    for (coupling_node, coupling_element) in pairs(coupling_dict)
        topo = vcat(coupling_node, topology[coupling_element])
        @views local_disp = displacements[topo, :]
        force_densities +=
            coupling_matrix[coupling_node, :, :] * local_disp * weight_coefficient
    end

    # Add weightening
    # for ipd in coupling_nodes
    #     force_internal[ipd*dof-1:ipd*dof] *= weight_coefficient
    # end
    #
    # for ife in coupling_FE_nodes
    #     force_internal[ife*dof-1:ife*dof] *= (1 - weight_coefficient)
    # end

    # Add weightening for mass vector same way

    # Add to the system coupling matrix
    #for (point_idx, point_val) in enumerate(coupling_nodes)
    #    topo_local = [point_val; topology[point_locations[point_idx], :]] # global numbering of coupled nodes [xpd x1fe x2fe x3fe x4fe]
    #
    #    # depends how vectors are constructed
    #end


    return datamanager

end




function compute_coupling_matrix(
    coordinates,
    topology,
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
    ### Right now only for linear elements
    x1 = coordinates[topology[element_number, 1], 1]
    y1 = coordinates[topology[element_number, 1], 2]
    x2 = coordinates[topology[element_number, 2], 1]
    y2 = coordinates[topology[element_number, 2], 2]
    x3 = coordinates[topology[element_number, 3], 1]
    y3 = coordinates[topology[element_number, 3], 2]
    x4 = coordinates[topology[element_number, 4], 1]
    y4 = coordinates[topology[element_number, 4], 2]
    ###
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
    N1p = 0.25 * Nxi[1] * Neta[1] #(1 - ksi) * (1 - eta)
    N2p = 0.25 * Nxi[2] * Neta[1] #(1 + ksi) * (1 - eta)
    N3p = 0.25 * Nxi[2] * Neta[2] #(1 + ksi) * (1 + eta)
    N4p = 0.25 * Nxi[1] * Neta[2] #(1 - ksi) * (1 + eta)

    I = ones(1, 1)
    Np = [N1p N2p N3p N4p]

    local_coupl_matrix = kappa * [I -Np; -Np' Np'*Np]

    return local_coupl_matrix
end


function find_point_in_elements(coordinates, topology, points_to_check)

    # Will store which element each point is in ->
    point_locations = Dict{Int64,Int64}()

    # Iterate through points to check
    @views for (point_index, point_coords) in
               enumerate(eachrow(coordinates[points_to_check, :]))
        # Check each element
        for (element_index, element_vertices) in enumerate(eachrow(topology))
            # Get vertex coordinates for this element
            poly_vertices = [coordinates[v, :] for v in element_vertices]

            # Check if point is in this element
            if is_point_in_polygon(point_coords, poly_vertices)
                point_locations[points_to_check[point_index]] = element_index
                break
            end
        end
    end

    return point_locations
end


function is_point_in_polygon(point, polygon_vertices)
    # Point-in-polygon test using ray casting method
    x, y = point
    n = length(polygon_vertices)
    inside = false

    j = n
    for i = 1:n
        @views xi, yi = polygon_vertices[i]
        @views xj, yj = polygon_vertices[j]

        intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)

        if intersect
            inside = !inside
        end

        j = i
    end

    return inside
end
#=

kappa = 1

weight_coefficient = 0.5 # just here, should be changed in future to a function




stiffPD = 0.01     # some factor just for simple case

# Time integration Solver

displacement = Array{Float64}(zeros(length(coordinates)))   # depends what structure has this vector in original code, here it is (x1,y1,x2,y2... )'

#displacement = fill(0.1, 10, dof)

for tt in 1:100

    force_internal_FE = stiffmatrFE * displacement[nodesFE*dof]
    force_internal_PD = force_dencities #here both x and Y coord
    force_internal = vcat(force_internal_FE, force_internal_PD)
    # calculation of damage

    # Add weightening
    for ipd in coupling_nodes
        force_internal[ipd*dof-1:ipd*dof] *= weight_coefficient
    end

    for ife in coupling_FE_nodes
        force_internal[ife*dof-1:ife*dof] *= weight_coefficient
    end

    # Add weightening for mass vector same way

    # Add to the system coupling matrix
    for (point_idx, point_val) in enumerate(coupling_nodes)
        topo_local = [point_val; topology[point_locations[point_idx], :]] # global numbering of coupled nodes [xpd x1fe x2fe x3fe x4fe]
        ...
        # depends how vectors are constructed
    end
end


=#

end
