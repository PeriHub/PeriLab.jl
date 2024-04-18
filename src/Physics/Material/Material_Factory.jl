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
include("./Material_Models/Correspondence.jl")
using Reexport
@reexport using .Correspondence
export init_material_model
export compute_forces
export determine_isotropic_parameter
export distribute_force_densities
export init_material_model_fields


"""
    init_material_model_fields(datamanager::Module)

Initialize material model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_material_model_fields(datamanager::Module)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Forces", Float64, dof) #-> only if it is an output
    # tbd later in the compute class
    datamanager.create_constant_node_field("External Forces", Float64, dof)
    datamanager.create_node_field("Force Densities", Float64, dof)
    datamanager.create_constant_node_field("External Force Densities", Float64, dof)
    datamanager.create_constant_node_field("Acceleration", Float64, dof)
    datamanager.create_node_field("Velocity", Float64, dof)
    datamanager.create_constant_bond_field("Bond Forces", Float64, dof)
    datamanager.set_synch("Bond Forces", false, false)
    datamanager.set_synch("Force Densities", true, false)
    datamanager.set_synch("Velocity", false, true)
    datamanager.set_synch("Displacements", false, true)
    datamanager.set_synch("Acceleration", false, true)
    datamanager.set_synch("Deformed Coordinates", false, true)

    return datamanager
end


"""
    init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the material model.

# Arguments
- `datamanager::Data_manager`: Datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `block::Int64`: Block.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)
    model_param = datamanager.get_properties(block, "Material Model")
    if !haskey(model_param, "Material Model")
        @error "Block " * string(block) * " has no material model defined."
        return nothing
    end

    if occursin("Correspondence", model_param["Material Model"])
        datamanager.set_model_module("Correspondence", Correspondence)
        datamanager = Correspondence.init_material_model(datamanager, nodes, model_param)
        return datamanager
    end

    material_models = split(model_param["Material Model"], "+")
    material_models = map(r -> strip(r), material_models)

    for material_model in material_models
        mod = Set_modules.create_module_specifics(material_model, module_list, "material_name")
        if isnothing(mod)
            @error "No material of name " * material_model * " exists."
        end
        datamanager.set_model_module(material_model, mod)
        datamanager = mod.init_material_model(datamanager, nodes, model_param)
    end
    #TODO in extra function
    nlist = datamanager.get_nlist()
    nlist_filtered_ids = datamanager.get_filtered_nlist()
    if !isnothing(nlist_filtered_ids)
        bond_norm = datamanager.get_field("Bond Norm")
        bond_geometry = datamanager.get_field("Bond Geometry")
        for iID in nodes
            if length(nlist_filtered_ids[iID]) != 0
                for neighborID in nlist_filtered_ids[iID]
                    bond_norm[iID][neighborID, :] .*= sign(dot((bond_geometry[iID][neighborID, :]), bond_norm[iID][neighborID, :]))
                end
            end
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
    material_models = split(model_param["Material Model"], "+")
    if occursin("Correspondence", model_param["Material Model"])
        mod = datamanager.get_model_module("Correspondence")
        datamanager = mod.compute_forces(datamanager, nodes, model_param, time, dt, to)
        return datamanager
    end
    for material_model in material_models
        mod = datamanager.get_model_module(material_model)
        datamanager = mod.compute_forces(datamanager, nodes, model_param, time, dt, to)
    end
    return datamanager
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
    nlist_filtered_ids = datamanager.get_filtered_nlist()
    bond_force = datamanager.get_field("Bond Forces")
    force_densities = datamanager.get_field("Force Densities", "NP1")
    volume = datamanager.get_field("Volume")
    bond_damage = datamanager.get_bond_damage("NP1")

    if !isnothing(nlist_filtered_ids)
        bond_norm = datamanager.get_field("Bond Norm")
        displacements = datamanager.get_field("Displacements", "NP1")
        force_densities = distribute_forces(nodes, nlist, nlist_filtered_ids, bond_force, volume, bond_damage, displacements, bond_norm, force_densities)
    else
        force_densities = distribute_forces(nodes, nlist, bond_force, volume, bond_damage, force_densities)
    end
    return datamanager
end

"""
    calculate_von_mises_stress(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Calculate the von Mises stress.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function calculate_von_mises_stress(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    dof = datamanager.get_dof()
    stress_NP1 = datamanager.get_field("Cauchy Stress", "NP1")
    von_Mises_stress = datamanager.get_field("von Mises Stress", "NP1")
    for iID in nodes
        von_Mises_stress[iID], spherical_stress_NP1, deviatoric_stress_NP1 = get_von_mises_stress(von_Mises_stress[iID], dof, stress_NP1[iID, :, :])
    end
    return datamanager
end

end
