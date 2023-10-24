# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_solver_name(params::Dict)
    if check_element(params["Solver"], "Verlet")
        return "Verlet"
    end
    return ""
end

function get_initial_time(params::Dict)

    if check_element(params["Solver"], "Initial Time")
        return params["Solver"]["Initial Time"]
    end

    @error "No initial time defined"
end

function get_final_time(params::Dict)

    if check_element(params["Solver"], "Final Time")
        return params["Solver"]["Final Time"]
    end
    @error "No final time defined"
end

function get_safety_factor(params::Dict)
    if check_element(params["Solver"]["Verlet"], "Safety Factor")
        return params["Solver"]["Verlet"]["Safety Factor"]
    end
    return 1.0
end

function get_fixed_dt(params::Dict)
    if check_element(params["Solver"]["Verlet"], "Fixed dt")
        return params["Solver"]["Verlet"]["Fixed dt"]
    end
    return true
end

function get_numerical_damping(params::Dict)
    if check_element(params["Solver"], "Numerical Damping")
        return params["Solver"]["Numerical Damping"]
    end
    return Float64(0.0)
end

function get_solver_options(params::Dict)
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