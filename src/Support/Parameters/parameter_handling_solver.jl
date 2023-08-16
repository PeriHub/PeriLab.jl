function get_solver_name(params)
    check = check_element(params["Solver"], "Verlet")
    if check_element(params["Solver"], "Verlet")

        return params["Solver"]["Verlet"]

    end
    return params["Discretization"]["Input Mesh File"]
end

function get_initial_time(params)
    check = check_element(params["Solver"], "Initial Time")
    if !check
        @error "No initial time defined"
    end
    return params["Discretization"]["Initial Time"]
end

function get_final_time(params)
    check = check_element(params["Solver"], "Final Time")
    if !check
        @error "No final time defined"
    end
    return params["Discretization"]["Final Time"]
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
    return Dict("Additive Models" => additive, "Damage Models" => damage, "Material Models" => mechanical, "Thermal Models" => thermal)
end