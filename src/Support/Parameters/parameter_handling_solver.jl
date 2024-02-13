# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export get_initial_time
export get_solver_options
export get_solver_name
export get_safety_factor
export get_final_time
export get_fixed_dt
export get_nsteps
export get_max_damage
export get_numerical_damping
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
    elseif haskey(params["Solver"], "External")
        return "External"
    end
    @error "Wrong or missing solvername. Verlet and External are the options."
    return nothing
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
        return Float64(params["Solver"]["Final Time"])
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
    return get(params["Solver"][get_solver_name(params)], "Safety Factor", 1.0)
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
    return get(params["Solver"][get_solver_name(params)], "Fixed dt", -1.0)
end


"""
    get_nsteps(params::Dict)

Get the fixed time step

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `nsteps::Int64`: The fixed time step
"""
function get_nsteps(params::Dict)
    return get(params["Solver"][get_solver_name(params)], "Number of Steps", 1)
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
    return get(params["Solver"], "Numerical Damping", Float64(0.0))
end

"""
get_max_damage(params::Dict)

Get the maximum damage.

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `write_after_damage::Bool`: The value
"""
function get_max_damage(params::Dict)
    return get(params["Solver"], "Maximum Damage", Inf64)
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
    additive::Bool = get(params["Solver"], "Additive Models", false)
    damage::Bool = get(params["Solver"], "Damage Models", false)
    mechanical::Bool = get(params["Solver"], "Material Models", true)
    thermal::Bool = get(params["Solver"], "Thermal Models", false)
    return Dict{String,Any}("Additive Models" => additive, "Damage Models" => damage, "Material Models" => mechanical, "Thermal Models" => thermal)
end