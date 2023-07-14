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