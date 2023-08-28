function get_solver_name(params)
    if check_element(params["Solver"], "Verlet")
        return "Verlet"
    end
    return ""
end

function get_initial_time(params)

    if check_element(params["Solver"], "Initial Time")
        return params["Solver"]["Initial Time"]
    end

    @error "No initial time defined"
end

function get_final_time(params)

    if check_element(params["Solver"], "Final Time")
        return params["Solver"]["Final Time"]
    end
    @error "No final time defined"
end

function get_safety_factor(params)
    if check_element(params["Solver"]["Verlet"], "Safety Factor")
        return params["Solver"]["Verlet"]["Safety Factor"]
    end
    return 1.0
end

function get_fixed_dt(params)
    if check_element(params["Solver"]["Verlet"], "Fixed dt")
        return params["Solver"]["Verlet"]["Fixed dt"]
    end
    return -1.0
end

function get_numerical_damping(params)
    if check_element(params["Solver"], "Numerical Damping")
        return params["Solver"]["Numerical Damping"]
    end
    return 0.0
end

function get_solver_options(params)
    additive::Bool = false
    damage::Bool = false
    mechanical::Bool = true
    thermal::Bool = false
    if check_element(params["Solver"], "Material Models")
        mechanical = params["Solver"]["Material Models"]
    end
    if check_element(params["Solver"], "Thermal Models")
        thermal = params["Solver"]["Thermal Models"]
    end
    if check_element(params["Solver"], "Additive Models")
        additive = params["Solver"]["Additive Models"]
    end
    if check_element(params["Solver"], "Damage Models")
        damage = params["Solver"]["Damage Models"]
    end
    return Dict{String,Any}("Additive Models" => additive, "Damage Models" => damage, "Material Models" => mechanical, "Thermal Models" => thermal)
end