# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export get_solver_steps
export get_solver_params
export get_initial_time
export get_model_options
export get_solver_name
export get_safety_factor
export get_final_time
export get_fixed_dt
export get_max_damage
export get_numerical_damping
export get_nsteps
export get_calculation_options

"""
    get_solver_steps(params::Dict)

Get the solver steps

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `solver_steps::List`: The solver steps
"""
function get_solver_steps(params::Dict)
    step_id = []
    for step_name in keys(params["Solver"])
        if typeof(params["Solver"][step_name]) == Dict{Any,Any} &&
           haskey(params["Solver"][step_name], "Step ID")
            append!(step_id, params["Solver"][step_name]["Step ID"])
        end
    end
    if length(step_id) == 0
        return [nothing]
    end
    return sort!(step_id)
end

"""
    get_solver_params(params::Dict)

Get the solver parameters

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `solver_params::Dict`: The solver parameters
"""
function get_solver_params(params::Dict, step_id)
    for step_name in keys(params["Solver"])
        if params["Solver"][step_name]["Step ID"] == step_id
            return params["Solver"][step_name]
        end
    end
    @error "Step ID $step_id not found"
    return nothing
end

"""
    get_solver_name(params::Dict)

Get the name of the solver

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `solver_name::String`: The name of the solver
"""
function get_solver_name(params::Dict)
    if haskey(params, "Verlet")
        return "Verlet"
    elseif haskey(params, "External")
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

    if haskey(params, "Initial Time")
        return Float64(params["Initial Time"])
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

    if haskey(params, "Final Time")
        return Float64(params["Final Time"])
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
    return Float64(get(params[get_solver_name(params)], "Safety Factor", 1.0))
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
    return Float64(get(params[get_solver_name(params)], "Fixed dt", -1.0))
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
    return get(params[get_solver_name(params)], "Number of Steps", 1)
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
    return Float64(get(params, "Numerical Damping", 0.0))
end

"""
get_max_damage(params::Dict)

Get the maximum damage.

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `max_damage::Float64`: The value
"""
function get_max_damage(params::Dict)
    return Float64(get(params, "Maximum Damage", Inf64))
end

"""
    get_model_options(params::Dict)

Get the solver options

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `solver_options::Dict`: The solver options
"""
function get_model_options(params::Dict)
    additive::Bool = get(params, "Additive Models", false)
    corrosion::Bool = get(params, "Corrosion Models", false)
    damage::Bool = get(params, "Damage Models", false)
    material::Bool = get(params, "Material Models", true)
    thermal::Bool = get(params, "Thermal Models", false)
    pre_calculation::Bool = get(params, "Pre Calculation Models", true)
    ###
    # TODO here, the order is predefined. Wrong place!
    ###
    return Vector{String}(
        filter(
            !isnothing,
            [
                if additive
                    "Additive"
                end
                if damage
                    "Damage"
                end
                if pre_calculation
                    "Pre_Calculation"
                end
                if thermal
                    "Thermal"
                end
                if corrosion
                    "Corrosion"
                end
                if material
                    "Material"
                end
            ],
        ),
    )

end

"""
    get_calculation_options(params::Dict)

Get the calculation options

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `solver_options::Dict`: The solver options
"""
function get_calculation_options(params::Dict)
    cauchy::Bool = get(params, "Calculate Cauchy", false)
    von_mises::Bool = get(params, "Calculate von Mises stress", false)
    strain::Bool = get(params, "Calculate Strain", false)
    return Dict{String,Any}(
        "Calculate Cauchy" => cauchy,
        "Calculate von Mises stress" => von_mises,
        "Calculate Strain" => strain,
    )
end
