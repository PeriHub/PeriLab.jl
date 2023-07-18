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
    mechanical::Bool = true
    additive::Bool = false
    thermal::Bool = false
    if check_element(params["Solver"], "Mechanical")
        mechanical = params["Solver"]["Mechanical"]
    end
    if check_element(params["Solver"], "Thermal")
        thermal = params["Solver"]["Thermal"]
    end
    if check_element(params["Solver"], "Additive")
        additive = params["Solver"]["Additive"]
    end
    return mechanical, thermal, additive
end