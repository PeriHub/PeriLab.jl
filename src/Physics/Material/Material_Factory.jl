# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Material
include("../../Core/Module_inclusion/set_Modules.jl")
include("material_basis.jl")
using .Set_modules
using TimerOutputs

global module_list = Set_modules.find_module_files(@__DIR__, "material_name")
Set_modules.include_files(module_list)

export init_material_model
export compute_forces
export determine_isotropic_parameter
export distribute_force_densities

"""
    init_material_model(datamanager::Module, block::Int64)

Initializes the material model.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `block::Int64`: Block.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module, block::Int64)
    model_param = datamanager.get_properties(block, "Material Model")
    specifics = Dict{String,String}("Call Function" => "init_material_model", "Name" => "material_name")
    if !haskey(model_param, "Material Model")
        @error "Block " * string(block) * " has no material model defined."
    end
    material_models = split(model_param["Material Model"], "+")
    material_models = map(r -> strip(r), material_models)

    for material_model in material_models
        datamanager = Set_modules.create_module_specifics(material_model, module_list, specifics, (datamanager,))
        if isnothing(datamanager)
            @error "No material of name " * material_model * " exists."
        end

    end
    return datamanager
end

"""
    compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64, to::TimerOutput)

Compute the forces.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `model_param::Dict`: The material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64, to::TimerOutput)
    specifics = Dict{String,String}("Call Function" => "compute_forces", "Name" => "material_name")
    material_models = split(model_param["Material Model"], "+")

    for material_model in material_models
        datamanager = Set_modules.create_module_specifics(material_model, module_list, specifics, (datamanager, nodes, model_param, time, dt, to))
        if isnothing(datamanager)
            @error "No material of name " * material_model * " exists."
        end

    end
    return datamanager
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
    distribute_force_densities(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Distribute the force densities.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function distribute_force_densities(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    nlist = datamanager.get_nlist()
    bond_force = datamanager.get_field("Bond Forces")
    force_densities = datamanager.get_field("Force Densities", "NP1")
    volume = datamanager.get_field("Volume")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    force_densities = distribute_forces(nodes, nlist, bond_force, volume, bond_damage, force_densities)
    return datamanager
end

end
