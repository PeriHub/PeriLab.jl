function get_number_of_materials(params)
    return length(params["Materials"])
end

function get_materials(params)
    return params["Materials"]
end
function get_material_names(params)
    return keys(params["Materials"])
end

function get_material(params, name)
    return params["Materials"][name]
end
function get_material_types(params)
    if check_element(params["Materials"], "Material Types")
        return params["Materials"]["Material Types"]
    elseif check_element(params, "Material Types")
        return check_element(params, "Material Types")
    else
        return Dict()
    end
end