# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Material
include("../../Core/Module_inclusion/set_Modules.jl")
include("Material_Models/Ordinary/Ordinary.jl")
using ...Material_Basis:
                         get_all_elastic_moduli,
                         distribute_forces!,
                         check_symmetry,
                         get_all_elastic_moduli,
                         init_local_damping_due_to_damage,
                         local_damping_due_to_damage
using LinearAlgebra: dot
# using .Set_modules
using TimerOutputs
using StaticArrays

global module_list = find_module_files(@__DIR__, "material_name")
include_files(module_list)
export init_model
export compute_model
export determine_isotropic_parameter
export distribute_force_densities
export init_fields
export fields_for_local_synchronization
export compute_local_damping
export init_local_damping

function compute_local_damping(datamanager::Module, nodes, params, dt)
    return local_damping_due_to_damage(datamanager, nodes, params, dt)
end
function init_local_damping(datamanager, nodes, material_parameter, damage_parameter)
    return init_local_damping_due_to_damage(datamanager,
                                            nodes,
                                            material_parameter,
                                            damage_parameter)
end

"""
    init_fields(datamanager::Module)

Initialize material model fields

# Arguments
- `datamanager::Data_Manager`: Datamanager.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function init_fields(datamanager::Module)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Forces", Float64, dof) #-> only if it is an output
    # tbd later in the compute class
    datamanager.create_constant_node_field("External Forces", Float64, dof)
    datamanager.create_node_field("Force Densities", Float64, dof)
    datamanager.create_constant_node_field("External Force Densities", Float64, dof)
    datamanager.create_constant_node_field("Acceleration", Float64, dof)
    datamanager.create_node_field("Velocity", Float64, dof)
    datamanager.create_constant_bond_field("Bond Forces", Float64, dof)
    datamanager.create_constant_bond_field("Temporary Bond Field", Float64, 1)
    deformed_coorN,
    deformed_coorNP1 = datamanager.create_node_field("Deformed Coordinates",
                                                     Float64, dof)
    deformed_coorN = copy(datamanager.get_field("Coordinates"))
    deformed_coorNP1 = copy(datamanager.get_field("Coordinates"))
    datamanager.create_node_field("Displacements", Float64, dof)
    datamanager.create_bond_field("Deformed Bond Geometry", Float64, dof)
    datamanager.create_bond_field("Deformed Bond Length", Float64, 1)
    # datamanager.set_synch("Bond Forces", false, true)
    datamanager.set_synch("Force Densities", true, false)
    datamanager.set_synch("Velocity", false, true)
    datamanager.set_synch("Displacements", false, true)
    datamanager.set_synch("Acceleration", false, true)
    datamanager.set_synch("Deformed Coordinates", false, true)

    return datamanager
end

"""
    init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the material model.

# Arguments
- `datamanager::Data_Manager`: Datamanager
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function init_model(datamanager::Module, nodes::AbstractVector{Int64},
                    block::Int64)
    model_param = datamanager.get_properties(block, "Material Model")::Dict{String,Any}
    if !haskey(model_param, "Material Model")
        @error "Block " * string(block) * " has no material model defined."
        return nothing
    end

    if occursin("Correspondence", model_param["Material Model"])
        datamanager.set_model_module("Correspondence", Correspondence)
        return Correspondence.init_model(datamanager, nodes, block, model_param)
    end

    material_models = split(model_param["Material Model"], "+")
    material_models = map(r -> strip(r), material_models)

    for material_model in material_models
        mod = create_module_specifics(material_model,
                                      module_list,
                                      "material_name")
        datamanager.set_analysis_model("Material Model", block, material_model)
        if isnothing(mod)
            @error "No material of name " * material_model * " exists."
        end
        datamanager.set_model_module(material_model, mod)
        datamanager = mod.init_model(datamanager, nodes, model_param)
    end
    #TODO in extra function
    # nlist = datamanager.get_nlist()
    nlist_filtered_ids = datamanager.get_filtered_nlist()
    if !isnothing(nlist_filtered_ids)
        bond_norm = datamanager.get_field("Bond Norm")
        bond_geometry = datamanager.get_field("Bond Geometry")
        for iID in nodes
            if length(nlist_filtered_ids[iID]) != 0
                for neighborID in nlist_filtered_ids[iID]
                    bond_norm[iID][neighborID] .*= sign(dot((bond_geometry[iID][neighborID]),
                                                            bond_norm[iID][neighborID]))
                end
            end
        end
    end

    return datamanager
end

"""
    fields_for_local_synchronization(datamanager, model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `datamanager::Module`: datamanager.
- `model::String`: Model class.
- `block::Int64`: block ID
# Returns
- `datamanager::Module`: Datamanager.
"""
function fields_for_local_synchronization(datamanager, model, block)
    model_param = datamanager.get_properties(block, "Material Model")
    if occursin("Correspondence", model_param["Material Model"])
        return Correspondence.fields_for_local_synchronization(datamanager, model, block,
                                                               model_param)
    end

    for material_model in datamanager.get_analysis_model("Material Model", block)
        mod = datamanager.get_model_module(material_model)
        mod.fields_for_local_synchronization(datamanager, model)
    end
    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, model_param::Dict{String,Any}, block::Int64, time::Float64, dt::Float64,to::TimerOutput,)

Computes the material models

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict{String,Any}`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       model_param::Dict{String,Any},
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    @timeit to "all" begin
        if occursin("Correspondence", model_param["Material Model"])
            @timeit to "corresponcence" begin
                datamanager = Correspondence.compute_model(datamanager, nodes, model_param,
                                                           block, time,
                                                           dt,
                                                           to)
                return datamanager
            end
        end

        for material_model in datamanager.get_analysis_model("Material Model", block)
            mod = datamanager.get_model_module(material_model)

            datamanager = mod.compute_model(datamanager, nodes, model_param, block, time,
                                            dt,
                                            to)
        end
        return datamanager
    end
end

"""
    determine_isotropic_parameter(datamanager::Module, prop::Dict)

Determine the isotropic parameter.

# Arguments
- `prop::Dict`: The material property.
# Returns
- `prop::Dict`: The material property.
"""
function determine_isotropic_parameter(datamanager::Module, prop::Dict)
    get_all_elastic_moduli(datamanager, prop)
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
    distribute_force_densities(datamanager::Module, nodes::AbstractVector{Int64})

Distribute the force densities.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: The nodes.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function distribute_force_densities(datamanager::Module,
                                    nodes::AbstractVector{Int64}, to::TimerOutput)
    @timeit to "load data" begin
        nlist = datamanager.get_nlist()
        nlist_filtered_ids = datamanager.get_filtered_nlist()
        bond_force = datamanager.get_field("Bond Forces")
        force_densities = datamanager.get_field("Force Densities", "NP1")
        volume = datamanager.get_field("Volume")
        bond_damage = datamanager.get_bond_damage("NP1")
    end
    if !isnothing(nlist_filtered_ids)
        bond_norm = datamanager.get_field("Bond Norm")
        displacements = datamanager.get_field("Displacements", "NP1")
        @timeit to "local dist" force_densities=distribute_forces!(force_densities,
                                                                   nodes,
                                                                   nlist,
                                                                   nlist_filtered_ids,
                                                                   bond_force,
                                                                   volume,
                                                                   bond_damage,
                                                                   displacements,
                                                                   bond_norm)
    else
        @timeit to "local dist" distribute_forces!(force_densities, nodes, nlist,
                                                   bond_force, volume, bond_damage)
    end
end
end
