# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Material

include("Material_Models/Ordinary/Ordinary.jl")

using TimerOutputs: @timeit
using ....Data_Manager
using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "material_name")
for mod in module_list
    include(mod["File"])
end

using ...Material_Basis:
                         get_all_elastic_moduli,
                         distribute_forces!,
                         check_symmetry,
                         get_all_elastic_moduli,
                         init_local_damping_due_to_damage,
                         local_damping_due_to_damage
using LinearAlgebra: dot
using StaticArrays

export init_model
export compute_model
export determine_isotropic_parameter
export distribute_force_densities
export init_fields
export fields_for_local_synchronization
export compute_local_damping
export init_local_damping

function compute_local_damping(nodes, params, dt)
    return local_damping_due_to_damage(nodes, params, dt)
end
function init_local_damping(nodes, material_parameter, damage_parameter)
    return init_local_damping_due_to_damage(nodes,
                                            material_parameter,
                                            damage_parameter)
end

"""
    init_fields()

Initialize material model fields
"""
function init_fields()
    dof = Data_Manager.get_dof()
    Data_Manager.create_node_field("Forces", Float64, dof) #-> only if it is an output
    # tbd later in the compute class
    Data_Manager.create_constant_node_field("External Forces", Float64, dof)
    Data_Manager.create_node_field("Force Densities", Float64, dof)
    Data_Manager.create_constant_node_field("External Force Densities", Float64, dof)
    Data_Manager.create_constant_node_field("Acceleration", Float64, dof)
    Data_Manager.create_node_field("Velocity", Float64, dof)
    Data_Manager.create_constant_bond_field("Bond Forces", Float64, dof)
    Data_Manager.create_constant_bond_field("Temporary Bond Field", Float64, 1)
    deformed_coorN,
    deformed_coorNP1 = Data_Manager.create_node_field("Deformed Coordinates",
                                                      Float64, dof)
    deformed_coorN = copy(Data_Manager.get_field("Coordinates"))
    deformed_coorNP1 = copy(Data_Manager.get_field("Coordinates"))
    Data_Manager.create_node_field("Displacements", Float64, dof)
    Data_Manager.create_bond_field("Deformed Bond Geometry", Float64, dof)
    Data_Manager.create_bond_field("Deformed Bond Length", Float64, 1)
    # Data_Manager.set_synch("Bond Forces", false, true)
    Data_Manager.set_synch("Force Densities", true, false)
    Data_Manager.set_synch("Velocity", false, true)
    Data_Manager.set_synch("Displacements", false, true)
    Data_Manager.set_synch("Acceleration", false, true)
    Data_Manager.set_synch("Deformed Coordinates", false, true)
end

"""
    init_model(nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the material model.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    model_param = Data_Manager.get_properties(block, "Material Model")::Dict{String,Any}
    if !haskey(model_param, "Material Model")
        @error "Block " * string(block) * " has no material model defined."
        return nothing
    end

    if occursin("Correspondence", model_param["Material Model"])
        Data_Manager.set_model_module("Correspondence", Correspondence)
        return Correspondence.init_model(nodes, block, model_param)
    end

    material_models = split(model_param["Material Model"], "+")
    material_models = map(r -> strip(r), material_models)

    for material_model in material_models
        mod = create_module_specifics(material_model,
                                      module_list,
                                      @__MODULE__,
                                      "material_name")
        Data_Manager.set_analysis_model("Material Model", block, material_model)
        if isnothing(mod)
            @error "No material of name " * material_model * " exists."
        end
        Data_Manager.set_model_module(material_model, mod)
        mod.init_model(nodes, model_param)
    end
    #TODO in extra function
    # nlist = Data_Manager.get_nlist()
    nlist_filtered_ids = Data_Manager.get_filtered_nlist()
    if !isnothing(nlist_filtered_ids)
        bond_norm = Data_Manager.get_field("Bond Norm")
        bond_geometry = Data_Manager.get_field("Bond Geometry")
        for iID in nodes
            if length(nlist_filtered_ids[iID]) != 0
                for neighborID in nlist_filtered_ids[iID]
                    bond_norm[iID][neighborID] .*= sign(dot((bond_geometry[iID][neighborID]),
                                                            bond_norm[iID][neighborID]))
                end
            end
        end
    end
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block id
"""
function fields_for_local_synchronization(model, block)
    model_param = Data_Manager.get_properties(block, "Material Model")
    if occursin("Correspondence", model_param["Material Model"])
        return Correspondence.fields_for_local_synchronization(model, block,
                                                               model_param)
    end

    for material_model in Data_Manager.get_analysis_model("Material Model", block)
        mod = Data_Manager.get_model_module(material_model)
        mod.fields_for_local_synchronization(model)
    end
end

"""
    compute_model(nodes::AbstractVector{Int64}, model_param::Dict{String,Any}, block::Int64, time::Float64, dt::Float64)

Computes the material models

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict{String,Any}`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64},
                       model_param::Dict{String,Any},
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    @timeit "all" begin
        if occursin("Correspondence", model_param["Material Model"])
            @timeit "corresponcence" begin
                Correspondence.compute_model(nodes, model_param,
                                             block, time,
                                             dt)
                return
            end
        end

        for material_model in Data_Manager.get_analysis_model("Material Model", block)
            mod = Data_Manager.get_model_module(material_model)

            mod.compute_model(nodes, model_param, block, time,
                              dt)
        end
    end
end

"""
    determine_isotropic_parameter(prop::Dict)

Determine the isotropic parameter.

# Arguments
- `prop::Dict`: The material property.
# Returns
- `prop::Dict`: The material property.
"""
function determine_isotropic_parameter(prop::Dict)
    get_all_elastic_moduli(prop)
end

"""
    check_material_symmetry(dof::Int64, prop::Dict)

Check the symmetry of the material.

# Arguments
- `dof::Int64`: The degree of freedom.
- `prop::Dict`: The material property.
# Returns
- `prop::Dict`: The material property.
"""
function check_material_symmetry(dof::Int64, prop::Dict)
    return check_symmetry(prop, dof)
end

"""
    distribute_force_densities(nodes::AbstractVector{Int64})

Distribute the force densities.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
"""
function distribute_force_densities(nodes::AbstractVector{Int64})
    @timeit "load data" begin
        nlist = Data_Manager.get_nlist()
        nlist_filtered_ids = Data_Manager.get_filtered_nlist()
        bond_force = Data_Manager.get_field("Bond Forces")
        force_densities = Data_Manager.get_field("Force Densities", "NP1")
        volume = Data_Manager.get_field("Volume")
        bond_damage = Data_Manager.get_bond_damage("NP1")
    end
    if !isnothing(nlist_filtered_ids)
        bond_norm = Data_Manager.get_field("Bond Norm")
        displacements = Data_Manager.get_field("Displacements", "NP1")
        @timeit "local dist" force_densities=distribute_forces!(force_densities,
                                                                nodes,
                                                                nlist,
                                                                nlist_filtered_ids,
                                                                bond_force,
                                                                volume,
                                                                bond_damage,
                                                                displacements,
                                                                bond_norm)
    else
        @timeit "local dist" distribute_forces!(force_densities, nodes, nlist,
                                                bond_force, volume, bond_damage)
    end
end
end
