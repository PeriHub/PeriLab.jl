# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module FEM
include("../Core/Module_inclusion/set_Modules.jl")
include("./FEM_routines.jl")
#include("../Models/Material/material_basis.jl")
#using .Material_Basis: voigt_to_matrix
# in future using set modules for material
# test case is correspondence material
include("./../Models/Material/Material_Models/Correspondence/Correspondence_Elastic.jl")
using .Correspondence_Elastic


using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "element_name")
Set_modules.include_files(module_list)

function init_FEM(complete_params::Dict, datamanager::Module)
    if !haskey(complete_params, "FEM")
        return nothing
    end
    params = complete_params["FEM"]
    if isnothing(valid_models(params))
        return nothing
    end

    datamanager.set_properties("FEM", params)
    if !haskey(complete_params["Models"]["Material Models"], params["Material Model"])
        @error "The FEM material model $(params["Material Model"]) is not defined"
        return nothing
    end
    datamanager.set_property(
        "FEM",
        "Material Model",
        complete_params["Models"]["Material Models"][params["Material Model"]],
    )

    dof = datamanager.get_dof()
    nelements = datamanager.get_num_elements()
    elements::Vector{Int64} = 1:nelements
    p = get_polynomial_degree(params, dof)
    coordinates = datamanager.get_field("Coordinates")
    if isnothing(p)
        return p
    end
    if dof != 2 && dof != 3
        @error "Degree of freedom = $dof is not supported, only 2 and 3."
        return nothing
    end
    num_int = get_number_of_integration_points(p, dof)
    N = datamanager.create_constant_free_size_field(
        "N Matrix",
        Float64,
        (prod(num_int), prod(p .+ 1) * dof, dof),
    )

    B_matrix = datamanager.create_constant_free_size_field(
        "B Matrix",
        Float64,
        (nelements, prod(num_int), prod(p .+ 1) * dof, 3 * dof - 3),
    )

    strainN, strainNP1 = datamanager.create_free_size_field(
        "Element Strain",
        Float64,
        (nelements, prod(num_int), 3 * dof - 3),
    )
    stressN, stressNP1 = datamanager.create_free_size_field(
        "Element Stress",
        Float64,
        (nelements, prod(num_int), 3 * dof - 3),
    )
    strain_increment = datamanager.create_constant_free_size_field(
        "Element Strain Increment",
        Float64,
        (nelements, prod(num_int), 3 * dof - 3),
    )

    specifics = Dict{String,String}(
        "Call Function" => "create_element_matrices",
        "Name" => "element_name",
    )
    # B_elem only temporary
    N[:], B_elem = create_element_matrices(
        dof,
        p,
        Set_modules.create_module_specifics(params["Element Type"], module_list, specifics),
    )
    if isnothing(N) || isnothing(B_matrix)
        return nothing
    end
    specifics =
        Dict{String,String}("Call Function" => "init_element", "Name" => "element_name")
    datamanager = Set_modules.create_module_specifics(
        params["Element Type"],
        module_list,
        specifics,
        (datamanager, elements, params, p),
    )

    elements = Vector{Int64}(1:nelements)
    topology = datamanager.get_field("FE Topology")
    if length(topology[1, :]) != prod(p .+ 1)
        @error "Size of topology and polynomial degree does not match."
        return nothing
    end
    jacobian = datamanager.create_constant_free_size_field(
        "Element Jacobi Matrix",
        Float64,
        (nelements, prod(num_int), dof, dof),
    )
    determinant_jacobian = datamanager.create_constant_free_size_field(
        "Element Jacobi Determinant",
        Float64,
        (nelements, prod(num_int)),
    )
    jacobian, determinant_jacobian = get_Jacobian(
        elements,
        dof,
        topology,
        coordinates,
        B_elem,
        jacobian,
        determinant_jacobian,
    )

    lumped_mass = datamanager.create_constant_node_field("Lumped Mass Matrix", Float64, 1)
    rho = datamanager.get_field("Density")
    lumped_mass =
        get_lumped_mass(elements, dof, topology, N, determinant_jacobian, rho, lumped_mass)
    B_matrix = create_B_matrix(elements, dof, B_elem, jacobian, B_matrix)

    datamanager = get_FEM_nodes(datamanager, topology)

    return datamanager

end

function valid_models(params::Dict)
    if haskey(params, "Additive Model")
        @warn "Additive models are not supported for FEM yet"
    end
    if haskey(params, "Damage Model")
        @warn "Damage models are not supported for FEM"

    end
    if haskey(params, "Thermal Model")
        @warn "Thermal models are not supported for FEM yet"
    end
    if !haskey(params, "Material Model")
        @error "No material model has been defined for FEM in the block."
        return nothing
    else
        # in future -> FE support -> check with set modules
        if !Correspondence_Elastic.fe_support()
            @error "No FEM support for " * params["Material Model"]
            return nothing
        end
    end
    return true
end

function eval(
    datamanager::Module,
    elements::Union{SubArray,Vector{Int64}},
    params::Dict,
    time::Float64,
    dt::Float64,
)
    return calculate_FEM(
        datamanager,
        elements,
        params,
        Correspondence_Elastic.compute_stresses,
        time,
        dt,
    )
end

function get_FEM_nodes(datamanager::Module, topology::Matrix{Int64})
    fem_nodes = datamanager.create_constant_node_field("FE Nodes", Bool, 1)
    for el_topo in eachrow(topology)
        fem_nodes[el_topo] .= true
    end
    return datamanager
end


end
