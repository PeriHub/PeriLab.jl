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
    if check_element(params["Solver"], "Material")
        mechanical = params["Solver"]["Material"]
    end
    if check_element(params["Solver"], "Thermal")
        thermal = params["Solver"]["Thermal"]
    end
    if check_element(params["Solver"], "Additive")
        additive = params["Solver"]["Additive"]
    end
    if check_element(params["Solver"], "Damage")
        damage = params["Solver"]["Damage"]
    end
    return Dict("Additive" => additive, "Damage" => damage, "Material" => mechanical, "Thermal" => thermal)
end