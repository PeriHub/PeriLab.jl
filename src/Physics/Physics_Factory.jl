module Physics
include("./Additive/Additive_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("../Support/Parameters/parameter_handling.jl")
import .Additive
import .Damage
import .Material
import .Thermal
export get_physics


function get_physics(params)
    Additive.get_additive(params)
    Damage.get_damage(params)
    Material.get_material(params)
    Thermal.get_thermal(params)
    return
end

function read_properties(params, datamanager)
    datamanager.init_property()
    blocks = datamanager.get_block_list()
    prop_keys = datamanager.init_property()
    for block in blocks
        get_block_model_definition(params, block, prop_keys, datamanager.set_properties)
    end
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

end