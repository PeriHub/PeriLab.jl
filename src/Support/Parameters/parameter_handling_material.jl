function get_number_of_materials(params)
    check = check_element(params, "Materials")
    if check
        return length(params["Materials"])
    else
        @error "No materials defined"
    end
    return 0
end

function get_materials(params)
    check = check_element(params, "Materials")
    if !check
        @error "No materials defined"
    end

    return params["Materials"]
end
