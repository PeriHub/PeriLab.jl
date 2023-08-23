function get_model_parameter(params, model, id)
    if check_element(params["Physics"], model) == false
        @error model * " is defined in blocks, but no " * model * " defintion block exists"
        return Dict()
    end
    if check_element(params[model], id)
        return params["Physics"][model][id]
    else
        @error model * " model with name " * id * " is defined in blocks,but missing in the " * model * " definition."
        return Dict()
    end
end

