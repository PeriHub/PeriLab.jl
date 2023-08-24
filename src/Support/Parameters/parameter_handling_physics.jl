function get_model_parameter(params, model, id)
    if check_element(params["Physics"], model * "s") == false
        @error model * " is defined in blocks, but no " * model * "s definition block exists"
        return Dict()
    end
    if check_element(params["Physics"][model*"s"], id)
        return params["Physics"][model*"s"][id]
    else
        @error model * " model with name " * id * " is defined in blocks, but missing in the " * model * "s definition."
        return Dict()
    end
end

