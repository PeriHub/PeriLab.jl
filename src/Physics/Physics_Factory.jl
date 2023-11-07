# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Physics
include("../Support/helpers.jl")
include("./Additive/Additive_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("./Pre_calculation/Pre_Calculation_Factory.jl")
include("../Support/Parameters/parameter_handling.jl")

using .Additive
using .Damage
using .Material
using .Pre_calculation
using .Thermal
using TimerOutputs
export compute_models
export init_models
export read_properties
export init_additive_model_fields
export init_damage_model_fields
export init_material_model_fields
export init_thermal_model_fields


function compute_models(datamanager::Module, block_nodes::Union{SubArray,Vector{Int64}}, block::Int64, dt::Float64, time::Float64, options::Dict, synchronise_field, to::TimerOutput)

    if options["Additive Models"]
        if datamanager.check_property(block, "Additive Model")
            @timeit to "compute_additive_model" datamanager = Additive.compute_additive(datamanager, block_nodes, datamanager.get_properties(block, "Additive Model"), time, dt)
        end
    end
    active = datamanager.get_field("Active")
    nodes = block_nodes[find_active(active[block_nodes])]
    update_list = datamanager.get_field("Update List")

    if options["Damage Models"]
        #tbd damage specific pre_calculation-> in damage template
        if datamanager.check_property(block, "Damage Model") && datamanager.check_property(block, "Material Model")
            datamanager = Damage.set_bond_damage(datamanager, nodes)
            @timeit to "compute_damage_pre_calculation" datamanager = compute_damage_pre_calculation(datamanager, options, nodes, block, synchronise_field, time, dt)
            @timeit to "compute_damage" datamanager = Damage.compute_damage(datamanager, nodes, datamanager.get_properties(block, "Damage Model"), block, time, dt)
            update_list = datamanager.get_field("Update List")
            update_nodes = view(nodes, find_active(update_list[nodes]))
            force_densities = datamanager.get_field("Force Densities", "NP1")
            force_densities[nodes] .= 0.0
        end
    end
    update_nodes = view(nodes, find_active(update_list[nodes]))
    @timeit to "pre_calculation" datamanager = Pre_calculation.compute(datamanager, update_nodes, datamanager.get_physics_options(), time, dt)

    if options["Thermal Models"]
        if datamanager.check_property(block, "Thermal Model")
            @timeit to "compute_thermal_model" datamanager = Thermal.compute_thermal_model(datamanager, nodes, datamanager.get_properties(block, "Thermal Model"), time, dt)
        end
    end

    if options["Material Models"]
        if datamanager.check_property(block, "Material Model")
            @timeit to "compute_bond_forces" datamanager = Material.compute_forces(datamanager, update_nodes, datamanager.get_properties(block, "Material Model"), time, dt)
            datamanager = Material.distribute_force_densities(datamanager, nodes)
        end
    end
    return datamanager

end
function compute_damage_pre_calculation(datamanager::Module, options::Dict, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64)

    datamanager = Pre_calculation.compute(datamanager, nodes, datamanager.get_physics_options(), time, dt)
    if options["Thermal Models"]
        datamanager = Thermal.compute_thermal_model(datamanager, nodes, datamanager.get_properties(block, "Thermal Model"), time, dt)
    end
    if options["Material Models"]
        datamanager = Material.compute_forces(datamanager, nodes, datamanager.get_properties(block, "Material Model"), time, dt)
    end
    datamanager = Material.distribute_force_densities(datamanager, nodes)
    synchronise_field(datamanager.get_comm(), datamanager.get_synch_fields(), datamanager.get_overlap_map(), datamanager.get_field, "Force DensitiesNP1", "download_from_cores")
    synchronise_field(datamanager.get_comm(), datamanager.get_synch_fields(), datamanager.get_overlap_map(), datamanager.get_field, "Force DensitiesNP1", "upload_to_cores")
    datamanager = Damage.compute_damage_pre_calculation(datamanager, nodes, block, datamanager.get_properties(block, "Damage Model"), synchronise_field, time, dt)
    update_list = datamanager.get_field("Update List")
    update_list[nodes] .= false
    return datamanager
end
function get_block_model_definition(params::Dict, blockID::Int64, prop_keys::Vector{String}, properties)
    # properties function from datamanager
    if check_element(params["Blocks"], "block_" * string(blockID))
        block = params["Blocks"]["block_"*string(blockID)]
        for model in prop_keys
            if check_element(block, model)
                properties(blockID, model, get_model_parameter(params::Dict, model, block[model]))
            end
        end
    end
    return properties
end

function init_material_model_fields(datamanager::Module)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Forces", Float64, dof) #-> only if it is an output
    # tbd later in the compute class
    datamanager.create_node_field("Forces", Float64, dof)
    datamanager.create_node_field("Force Densities", Float64, dof)
    datamanager.create_constant_node_field("Acceleration", Float64, dof)
    datamanager.create_node_field("Velocity", Float64, dof)
    datamanager.set_synch("Force Densities", true, false)
    datamanager.set_synch("Velocity", false, true)
    datamanager.set_synch("Displacements", false, true)
    datamanager.set_synch("Acceleration", false, true)
    datamanager.set_synch("Deformed Coordinates", false, true)

    return datamanager
end

function init_damage_model_fields(datamanager::Module)
    datamanager.create_node_field("Damage", Float64, 1)
    return datamanager
end


function init_models(params::Dict, datamanager::Module, allBlockNodes::Dict{Int64,Vector{Int64}}, solver_options::Dict)
    dof = datamanager.get_dof()
    defCoorN, defCoorNP1 = datamanager.create_node_field("Deformed Coordinates", Float64, dof)
    defCoorN[:] = copy(datamanager.get_field("Coordinates"))
    defCoorNP1[:] = copy(datamanager.get_field("Coordinates"))
    datamanager.create_node_field("Displacements", Float64, dof)
    if solver_options["Additive Models"]
        datamanager = Physics.init_additive_model_fields(datamanager)
        heatCapacity = datamanager.create_constant_node_field("Specific Heat Capacity", Float64, 1)
        heatCapacity = set_heatcapacity(params, allBlockNodes, heatCapacity) # includes the neighbors
    end
    if solver_options["Damage Models"]
        datamanager = Physics.init_damage_model_fields(datamanager)
        datamanager = Damage.init_interface_crit_values(datamanager, params)
    end
    if solver_options["Material Models"]
        datamanager = Physics.init_material_model_fields(datamanager)
    end
    if solver_options["Thermal Models"]
        datamanager = Physics.init_thermal_model_fields(datamanager)
        heatCapacity = datamanager.create_constant_node_field("Specific Heat Capacity", Float64, 1)
        heatCapacity = set_heatcapacity(params, allBlockNodes, heatCapacity) # includes the neighbors
    end


    return init_pre_calculation(datamanager, datamanager.get_physics_options())
end

function init_thermal_model_fields(datamanager::Module)
    datamanager.create_node_field("Temperature", Float64, 1)
    datamanager.create_node_field("Heat Flow", Float64, 1)
    datamanager.create_node_field("Specific Volume", Float64, 1)
    datamanager.create_constant_bond_field("Bond Heat Flow", Float64, 1)
    # if it is already initialized via mesh file no new field is created here
    datamanager.create_constant_node_field("Surface_nodes", Bool, 1)
    return datamanager
end

function init_additive_model_fields(datamanager::Module)
    if !("Activation_Time" in datamanager.get_all_field_keys())
        @error "'Activation_Time' is missing. Please define an 'Activation_Time' for each point in the mesh file."
    end
    # must be specified, because it might be that no temperature model has been defined
    datamanager.create_node_field("Temperature", Float64, 1)
    datamanager.create_node_field("Heat Flow", Float64, 1)
    bond_damageN = datamanager.get_field("Bond Damage", "N")
    bond_damageNP1 = datamanager.get_field("Bond Damage", "NP1")
    nnodes = datamanager.get_nnodes()
    if !("Active" in datamanager.get_all_field_keys())
        active = datamanager.create_constant_node_field("Active", Bool, 1)
        active .= false
        for iID in 1:nnodes
            bond_damageN[iID][:] .= 0
            bond_damageNP1[iID][:] .= 0
        end
    end
    return datamanager
end

function init_pre_calculation(datamanager::Module, options)
    return Pre_calculation.init_pre_calculation(datamanager, options)
end

function read_properties(params::Dict, datamanager::Module, material_model::Bool)
    datamanager.init_property()
    blocks = datamanager.get_block_list()
    prop_keys = datamanager.init_property()
    physics_options = datamanager.get_physics_options()
    datamanager.set_physics_options(get_physics_option(params::Dict, physics_options))

    for block in blocks
        get_block_model_definition(params::Dict, block, prop_keys, datamanager.set_properties)
    end
    if material_model
        dof = datamanager.get_dof()
        for block in blocks
            Material.check_material_symmetry(dof, datamanager.get_properties(block, "Material Model"))
            Material.determine_isotropic_parameter(datamanager.get_properties(block, "Material Model"))
        end
    end
    return datamanager
end

function set_heatcapacity(params::Dict, blockNodes::Dict, heatCapacity::SubArray)
    for block in eachindex(blockNodes)
        heatCapacity[blockNodes[block]] .= get_heatcapacity(params, block)
    end
    return heatCapacity
end

end