
include("../../Support/tools.jl")
include("../../Support/Parameters/parameter_handling.jl")
module Material
export get_material

function get_material(params)
    materials = get_materials(params)
    material_names = keys(material)
    #files = find_files_with_ending("./", "*.jl")
    # blocks
    # blocks -> mat
end


function material_type(data)
    ## function for specific pre-calculations
    if data.material.correspondence
        include("Correspondence.jl")
        data = shapeTensor(data)
        data = defGrad(data)
    end
    if data.material.correspondence
    end
    return data
end
end