module Physics
include("./Additive/Additive_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("../Support/Parameters/parameter_handling.jl")

using .Additive
using .Damage
using .Material
using .Thermal

export read_properties


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

function PD_type(params, datamanager)
    ## function for specific pre-calculations

    if datamanager.is_material_type("Bond-Based")
    end

    if datamanager.is_material_type("Ordinary")
    end

    if datamanager.is_material_type("Correspondence")
        include("Correspondence.jl")
        data = shapeTensor(data)
        data = defGrad(data)
    end
    if datamanager.is_material_type("Bond Associated")
    end

    return data
end

end