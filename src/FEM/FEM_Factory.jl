# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module FEM

using ...Data_Manager
using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "element_name")
for mod in module_list
    include(mod["File"])
end

include("./FEM_basis.jl")
using .FEM_Basis:
                  compute_FEM,
                  create_element_matrices,
                  get_Jacobian,
                  get_number_of_integration_points,
                  get_lumped_mass,
                  get_polynomial_degree,
                  create_B_matrix,
                  create_B_matrix
include("./Coupling/Coupling_Factory.jl")

# in future using set modules for material
# test case is correspondence material
using ...Helpers: fast_mul!, get_mapping
using ..Material_Basis: get_Hooke_matrix
# using .Correspondence_Elastic
using .Coupling
export init_FEM
export eval_FEM

function init_FEM(complete_params::Dict)
    if !haskey(complete_params, "FEM")
        return nothing
    end
    params = convert(Dict{String,Any}, complete_params["FEM"])
    if isnothing(valid_models(params))
        return nothing
    end

    Data_Manager.set_properties("FEM", params)
    if !haskey(complete_params["Models"]["Material Models"], params["Material Model"])
        @error "The FEM material model $(params["Material Model"]) is not defined"
        return nothing
    end
    @info "Initialize FEM"
    @warn "FEM Material models are set to all blocks. TODO must be changed in future."
    Data_Manager.set_property("FEM",
                              "Material Model",
                              complete_params["Models"]["Material Models"][params["Material Model"]])

    dof = Data_Manager.get_dof()
    nelements = Data_Manager.get_num_elements()
    elements::Vector{Int64} = 1:nelements
    p = get_polynomial_degree(params, dof)
    coordinates = Data_Manager.get_field("Coordinates")
    if isnothing(p)
        return p
    end
    if dof != 2 && dof != 3
        @error "Degree of freedom = $dof is not supported, only 2 and 3."
        return nothing
    end
    num_int = get_number_of_integration_points(p, dof)
    N = Data_Manager.create_constant_free_size_field("N Matrix",
                                                     Float64,
                                                     (prod(num_int), prod(p .+ 1) * dof,
                                                      dof))

    B_matrix = Data_Manager.create_constant_free_size_field("B Matrix",
                                                            Float64,
                                                            (nelements, prod(num_int),
                                                             prod(p .+ 1) * dof,
                                                             3 * dof - 3))

    strainN,
    strainNP1 = Data_Manager.create_free_size_field("Element Strain",
                                                    Float64,
                                                    (nelements, prod(num_int),
                                                     3 * dof - 3))
    stressN,
    stressNP1 = Data_Manager.create_free_size_field("Element Stress",
                                                    Float64,
                                                    (nelements, prod(num_int),
                                                     3 * dof - 3))
    strain_increment = Data_Manager.create_constant_free_size_field("Element Strain Increment",
                                                                    Float64,
                                                                    (nelements,
                                                                     prod(num_int),
                                                                     3 * dof - 3))

    specifics = Dict{String,String}("Call Function" => "create_element_matrices",
                                    "Name" => "element_name")
    # B_elem only temporary
    N[:],
    B_elem = create_element_matrices(dof,
                                     p,
                                     create_module_specifics(params["Element Type"],
                                                             module_list,
                                                             @__MODULE__,
                                                             specifics))
    if isnothing(N) || isnothing(B_matrix)
        return nothing
    end
    specifics = Dict{String,String}("Call Function" => "init_element",
                                    "Name" => "element_name")
    create_module_specifics(params["Element Type"],
                            module_list,
                            @__MODULE__,
                            specifics,
                            (elements, params, p))

    elements = Vector{Int64}(1:nelements)
    topology = Data_Manager.get_field("FE Topology")
    if length(topology[1, :]) != prod(p .+ 1)
        @error "Size of topology and polynomial degree does not match."
        return nothing
    end
    jacobian = Data_Manager.create_constant_free_size_field("Element Jacobi Matrix",
                                                            Float64,
                                                            (nelements, prod(num_int), dof,
                                                             dof))
    determinant_jacobian = Data_Manager.create_constant_free_size_field("Element Jacobi Determinant",
                                                                        Float64,
                                                                        (nelements,
                                                                         prod(num_int)))
    jacobian,
    determinant_jacobian = get_Jacobian(elements,
                                        dof,
                                        topology,
                                        coordinates,
                                        B_elem,
                                        jacobian,
                                        determinant_jacobian)

    lumped_mass = Data_Manager.create_constant_node_scalar_field("Lumped Mass Matrix",
                                                                 Float64)
    rho = Data_Manager.get_field("Density")
    lumped_mass = get_lumped_mass(elements, dof, topology, N, determinant_jacobian, rho,
                                  lumped_mass)
    B_matrix = create_B_matrix(elements, dof, B_elem, jacobian, B_matrix)

    get_FEM_nodes(topology)
    @info "End FEM init"
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
        # else
        #     # in future -> FE support -> check with set modules
        #     if !Correspondence_Elastic.fe_support()
        #         @error "No FEM support for " * params["Material Model"]
        #         return nothing
        #     end
    end
    return true
end

function compute_stresses(dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::Vector{Float64},
                          stress_N::Vector{Float64},
                          stress_NP1::Vector{Float64})
    hookeMatrix = get_Hooke_matrix(material_parameter,
                                   material_parameter["Symmetry"],
                                   dof)

    return hookeMatrix * strain_increment + stress_N
end

function eval_FEM(elements::AbstractVector{Int64},
                  params::Dict{String,Any},
                  time::Float64,
                  dt::Float64)
    return compute_FEM(elements,
                       params,
                       compute_stresses,
                       time,
                       dt)
end

function get_FEM_nodes(topology::Matrix{Int64})
    fem_nodes = Data_Manager.create_constant_node_scalar_field("FE Nodes", Bool)
    for el_topo in eachrow(topology)
        fem_nodes[el_topo] .= true
    end
end

"""
    force_densities(nodes)

Computes the force densities from the FEM nodes.

# Arguments
- `nodes::Vector{Int64}`: FEM nodes.
# Returns

"""
function force_densities(nodes::AbstractVector{Int64})
    volume = Data_Manager.get_field("Volume")
    forces = Data_Manager.get_field("Forces", "NP1")
    force_densities = Data_Manager.get_field("Force Densities", "NP1")
    #forces[nodes, :] ./= volume[nodes] # unclear why its needed here
    force_densities[nodes, :] = forces[nodes, :] ./ volume[nodes]
end

end
