
include("../../Support/tools.jl")
include("../../Support/Parameters/parameter_handling.jl")
module Material
export get_material


function material_types(params, datamanager)
    material_types = get_material_types(params)
    for material_type in material_types
        datamanager.set_material_type(material_type, material_types[material_type])
    end
    return datamanager.get_material_type()
end

function get_material(params)
    materials = get_materials(params)
    material_names = keys(material)
    #files = find_files_with_ending("./", "*.jl")
    # blocks
    # blocks -> mat
end


function material_type(params, datamanager)
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