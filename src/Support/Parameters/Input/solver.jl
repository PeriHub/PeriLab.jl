# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    SolverOptions

Options section of one solver. A solver module declares
`@params struct MyParams <: SolverOptions` (explicit solvers by convention with `safety_factor`,
`fixed_dt` (-1 if not given) and `numerical_damping`) and registers it with
`register_solver(name, MyParams)` in its `__init__()`; the input deck names the
solver by that key inside `Solver` or a `Multistep Solver` step.
"""
abstract type SolverOptions end

"""
    FallbackParams

The solver that takes over when a step of the main solver fails; it recomputes
the failed step and runs to the final time of the solver step.
"""
@params "Solver that takes over when a step fails; it recomputes the failed step from the last converged state and runs to the end of the solver step" struct FallbackParams
    "The solver that takes over a failed step"
    options::SolverOptions = one_of(:solver)
end

"Used for `Solver` and for every step of `Multistep Solver`."
@params "Time range, solver and the model categories evaluated in a solver step" struct SolverParams
    "Number of the step in a Multistep Solver (required there); steps run in ascending order and boundary conditions refer to it"
    step_id::Union{Nothing,Int64} = opt("Step ID"; default = nothing)
    "Start time of the step (required in Solver); a step of a Multistep Solver starts at the current time if this is not given or earlier"
    initial_time::Union{Nothing,Float64} = opt("Initial Time"; default = nothing,
                                               quantity = :time)
    "End time of the step (required in Solver); a step of a Multistep Solver may give Additional Time instead"
    final_time::Union{Nothing,Float64} = opt("Final Time"; default = nothing, quantity = :time)
    "Duration of the step, added to the current time; give this or Final Time, which takes precedence"
    additional_time::Union{Nothing,Float64} = opt("Additional Time"; default = nothing,
                                                  quantity = :time)
    "The solver of this step"
    options::SolverOptions = one_of(:solver)
    "Number of time steps; ignored when Fixed dt is set. 1 if not given"
    number_of_steps::Union{Nothing,Int64} = opt("Number of Steps"; default = nothing, min = 1,
                                                when = Dict("\$solver" => ["Static", "Newmark",
                                                "Linear Static Matrix Based",
                                                "Verlet Matrix Based"]))
    "Evaluates the material models of the blocks"
    material_models::Bool = opt("Material Models"; default = true)
    "Evaluates the damage models of the blocks"
    damage_models::Bool = opt("Damage Models"; default = false)
    "Evaluates the thermal models of the blocks"
    thermal_models::Bool = opt("Thermal Models"; default = false)
    "Evaluates the pre-calculations of the blocks"
    pre_calculation_models::Bool = opt("Pre Calculation Models"; default = true)
    "Evaluates the additive models of the blocks"
    additive_models::Bool = opt("Additive Models"; default = false)
    "Evaluates the degradation models of the blocks"
    degradation_models::Bool = opt("Degradation Models"; default = false)
    "Stops the step once the maximum damage of any node exceeds this value"
    maximum_damage::Float64 = opt("Maximum Damage"; default = Inf)
    "Computes the Cauchy stress of non-correspondence materials (correspondence materials always do), e.g. for output"
    calculate_cauchy::Bool = opt("Calculate Cauchy"; default = false)
    "Computes the von Mises stress, e.g. for output"
    calculate_von_mises_stress::Bool = opt("Calculate von Mises stress"; default = false)
    "Computes the strain of non-correspondence materials from the stress, e.g. for output"
    calculate_strain::Bool = opt("Calculate Strain"; default = false)
    fallback::Union{Nothing,FallbackParams} = opt("Fallback"; default = nothing)
end

function ParameterSpec.required_groups(::Type{SolverParams})
    return [ParameterSpec.required_group(["Final Time", "Additional Time"])]
end

"Registered name of the solver whose options are `o`."
function solver_name(o::SolverOptions)
    for (name, T) in ParameterSpec.registered_models(:solver)
        o isa T && return name
    end
    throw(ArgumentError("$(typeof(o)) is not a registered solver"))
end
solver_name(s::SolverParams) = solver_name(s.options)

"""
    start_time(s, current_time)

Start time of the solver step: `Initial Time`, or the current time if that is
later or `Initial Time` is not given (later steps of a multistep run).
"""
function start_time(s::SolverParams, current_time::Float64)
    if s.initial_time !== nothing
        return max(s.initial_time, current_time)
    end
    current_time != 0.0 && return current_time
    @abort "No initial time defined"
end

"""
    end_time(s, current_time)

End time of the solver step: `Final Time`, or `Additional Time` after the
current time.
"""
function end_time(s::SolverParams, current_time::Float64)
    s.final_time !== nothing && return s.final_time
    s.additional_time !== nothing && return current_time + s.additional_time
    @abort "No final time defined"
end

"Active model categories, in the order the models are evaluated."
function model_options(s::SolverParams)
    return [name
            for (name, used) in (("Additive", s.additive_models), ("Damage", s.damage_models),
                                 ("Pre_Calculation", s.pre_calculation_models),
                                 ("Thermal", s.thermal_models),
                                 ("Degradation", s.degradation_models),
                                 ("Material", s.material_models)) if used]
end

"Which stress and strain fields the models compute."
struct CalculationOptions
    cauchy::Bool
    von_mises::Bool
    strain::Bool
end

function calculation_options(s::SolverParams)
    return CalculationOptions(s.calculate_cauchy, s.calculate_von_mises_stress,
                              s.calculate_strain)
end

"Merges the options of several solver steps: a field is computed if any step asks for it."
function Base.:|(a::CalculationOptions, b::CalculationOptions)
    return CalculationOptions(a.cauchy | b.cauchy, a.von_mises | b.von_mises,
                              a.strain | b.strain)
end
