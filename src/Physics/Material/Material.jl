function Material(data)
    if data.material.correspondence
        include("Correspondence.jl")
        data = shapeTensor(data)
    end
    return data
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