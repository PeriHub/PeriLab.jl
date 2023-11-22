# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    get_solver_name(params::Dict)

Get the name of the solver

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `solver_name::String`: The name of the solver
"""
function get_solver_name(params::Dict)
    if haskey(params["Solver"], "Verlet")
        return "Verlet"
    end
    @error "Wrong or missing solvername. Verlet is the only option yet."
    return
end

"""
    get_initial_time(params::Dict)

Get the initial time

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `initial_time::Float64`: The initial time
"""
function get_initial_time(params::Dict)

    if haskey(params["Solver"], "Initial Time")
        return params["Solver"]["Initial Time"]
    end

    @error "No initial time defined"
end

"""
    get_final_time(params::Dict)

Get the final time

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `final_time::Float64`: The final time
"""
function get_final_time(params::Dict)

    if haskey(params["Solver"], "Final Time")
        return params["Solver"]["Final Time"]
    end
    @error "No final time defined"
end

"""
    get_safety_factor(params::Dict)

Get the safety factor

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `safety_factor::Float64`: The safety factor
"""
function get_safety_factor(params::Dict)
    if haskey(params["Solver"]["Verlet"], "Safety Factor")
        return params["Solver"]["Verlet"]["Safety Factor"]
    end
    return 1.0
end

"""
    get_fixed_dt(params::Dict)

Get the fixed time step

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `fixed_dt::Float64`: The fixed time step
"""
function get_fixed_dt(params::Dict)
    if haskey(params["Solver"]["Verlet"], "Fixed dt")
        return params["Solver"]["Verlet"]["Fixed dt"]
    end
    return true
end

"""
    get_numerical_damping(params::Dict)

Get the numerical damping

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `numerical_damping::Float64`: The numerical damping
"""
function get_numerical_damping(params::Dict)
    if haskey(params["Solver"], "Numerical Damping")
        return params["Solver"]["Numerical Damping"]
    end
    return Float64(0.0)
end

"""
    get_solver_options(params::Dict)
    
Get the solver options

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `solver_options::Dict`: The solver options
"""
function get_solver_options(params::Dict)
    additive::Bool = false
    damage::Bool = false
    mechanical::Bool = true
    thermal::Bool = false
    if haskey(params["Solver"], "Material Models")
        mechanical = params["Solver"]["Material Models"]
    end
    if haskey(params["Solver"], "Thermal Models")
        thermal = params["Solver"]["Thermal Models"]
    end
    if haskey(params["Solver"], "Additive Models")
        additive = params["Solver"]["Additive Models"]
    end
    if haskey(params["Solver"], "Damage Models")
        damage = params["Solver"]["Damage Models"]
    end
    return Dict{String,Any}("Additive Models" => additive, "Damage Models" => damage, "Material Models" => mechanical, "Thermal Models" => thermal)
end