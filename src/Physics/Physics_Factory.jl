module Physics
include("./Additive/Additive_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("../Support/geometry.jl")
#include("../Support/Parameters/parameter_handling.jl")
using .Additive
using .Damage
using .Material
using .Thermal
using .Geometry
export compute_models
export read_properties
export init_additive_model_fields
export init_damage_model_fields
export init_material_model_fields
export init_thermal_model_fields



function compute_models(datamanager, block, dt, time)

    datamanager = pre_calculation(datamanager, datamanager.get_physics_options())


    datamanager = Material.compute_forces(datamanager, datamanager.get_properties(block, "Material Model"), time, dt)

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
    datamanager.create_node_field("Forces", Float32, dof)
    defCoorN, defCoorNP1 = datamanager.create_node_field("Deformed Coordinates", Float32, dof)
    defCoorN = copy(datamanager.get_field("Coordinates"))
    defCoorNP1 = copy(datamanager.get_field("Coordinates"))
    datamanager.create_node_field("Displacements", Float32, dof)
    datamanager.create_constant_node_field("Acceleration", Float32, dof)
    datamanager.create_node_field("Velocity", Float32, dof)
    datamanager.set_synch("Force", true, false)
    datamanager.set_synch("Velocity", false, true)
    datamanager.set_synch("Displacements", false, true)
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

function pre_calculation(datamanager, options)
    ## function for specific pre-calculations
    # sub functions in material, damage, etc.
    dof = datamanager.get_dof()
    nnodes = datamanager.get_nnodes()
    nlist = datamanager.get_nlist()
    defCoor = datamanager.get_field("Deformed Coordinates", "NP1")
    if options["Calculate Deformed Bond Geometry"]
        bond_defN, bond_defNP1 = datamanager.create_bond_field("Deformed Bond Geometry", Float32, dof + 1)

        bond_defNP1 = Geometry.bond_geometry(nnodes, dof, nlist, defCoor, bond_defNP1)
    end
    if options["Calculate Shape Tensor"]

    end
    if options["Calculate Deformation Gradient"]

    end
    if options["Calculate Bond Associated Shape Tensor"]

    end
    if options["Calculate Bond Associated Deformation Gradient"]

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