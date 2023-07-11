

module Material
export get_material

function get_material(params)

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