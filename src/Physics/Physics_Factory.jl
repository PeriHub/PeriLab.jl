# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Physics
include("./Additive/Additive_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("../Support/geometry.jl")
include("../Support/Parameters/parameter_handling.jl")
using .Additive
using .Damage
using .Material
using .Thermal
using .Geometry
using TimerOutputs
export compute_models
export read_properties
export init_additive_model_fields
export init_damage_model_fields
export init_material_model_fields
export init_thermal_model_fields

function init_models(datamanager)
    return init_pre_calculation(datamanager, datamanager.get_physics_options())
end
function compute_models(datamanager, nodes, block, dt, time, options, to)

    @timeit to "pre_calculation" datamanager = pre_calculation(datamanager, nodes, datamanager.get_physics_options())

    if options["Damage Models"]
        @timeit to "compute_damage" datamanager = Damage.compute_damage(datamanager, nodes, datamanager.get_properties(block, "Damage Model"), time, dt)
    end
    if options["Material Models"]
        @timeit to "compute_forces" datamanager = Material.compute_forces(datamanager, nodes, datamanager.get_properties(block, "Material Model"), time, dt)
    end

    if options["Thermal Models"]
        @warn "Thermal Models not included yet"
    end
    if options["Additive Models"]
        @warn "Additive Models not included yet"
    end
    return datamanager

end

function get_block_model_definition(params, blockID, prop_keys, properties)
    # properties function from datamanager
    if check_element(params["Blocks"], "block_" * string(blockID))
        block = params["Blocks"]["block_"*string(blockID)]
        for model in prop_keys
            if check_element(block, model)
                properties(blockID, model, get_model_parameter(params, model, block[model]))
            end
        end
    end
    return properties
end

function init_material_model_fields(datamanager)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Forces", Float32, dof) #-> only if it is an output
    # tbd later in the compute class
    datamanager.create_node_field("Forces", Float32, dof)
    datamanager.create_node_field("Force Densities", Float32, dof)
    defCoorN, defCoorNP1 = datamanager.create_node_field("Deformed Coordinates", Float32, dof)
    defCoorN[:] = copy(datamanager.get_field("Coordinates"))
    defCoorNP1[:] = copy(datamanager.get_field("Coordinates"))
    datamanager.create_node_field("Displacements", Float32, dof)
    datamanager.create_constant_node_field("Acceleration", Float32, dof)
    datamanager.create_node_field("Velocity", Float32, dof)
    datamanager.set_synch("Force Densities", true, false)
    datamanager.set_synch("Velocity", false, true)
    datamanager.set_synch("Displacements", false, true)
    datamanager.set_synch("Acceleration", false, true)
    datamanager.set_synch("Deformed Coordinates", false, true)
    return datamanager
end

function init_damage_model_fields(datamanager)
    datamanager.create_node_field("Damage", Float32, 1)
    return datamanager
end

function init_thermal_model_fields(datamanager)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Temperature", Float32, 1)
    datamanager.create_node_field("Heat Flow", Float32, dof)
    return datamanager
end
function init_additive_model_fields(datamanager)
    datamanager.create_constant_node_field("Activated", Bool, 1)
    return datamanager
end

function init_pre_calculation(datamanager, options)

    dof = datamanager.get_dof()

    if options["Deformed Bond Geometry"]
        bond_defN, bond_defNP1 = datamanager.create_bond_field("Deformed Bond Geometry", Float32, dof + 1)
    end
    if options["Shape Tensor"]

    end
    if options["Deformation Gradient"]

    end
    if options["Bond Associated Shape Tensor"]

    end
    if options["Bond Associated Deformation Gradient"]

    end

    return datamanager
end

function pre_calculation(datamanager, nodes, options)
    ## function for specific pre-calculations
    # sub functions in material, damage, etc.
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    defCoor = datamanager.get_field("Deformed Coordinates", "NP1")
    if options["Deformed Bond Geometry"]
        bond_defNP1 = datamanager.get_field("Deformed Bond Geometry", "NP1")
        bond_defNP1 = Geometry.bond_geometry(nodes, dof, nlist, defCoor, bond_defNP1)
    end

    if options["Shape Tensor"]

    end
    if options["Deformation Gradient"]

    end
    if options["Bond Associated Shape Tensor"]

    end
    if options["Bond Associated Deformation Gradient"]

    end

    return datamanager
end

function read_properties(params, datamanager)
    datamanager.init_property()
    blocks = datamanager.get_block_list()
    prop_keys = datamanager.init_property()
    physics_options = datamanager.get_physics_options()
    datamanager.set_physics_options(get_physics_options(params, physics_options))
    for block in blocks
        get_block_model_definition(params, block, prop_keys, datamanager.set_properties)
    end
end

end