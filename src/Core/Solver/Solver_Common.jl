# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    Solver_Common

What every solver uses: the run of a solver step (`SolverRun`), what the run
writes to (`SolverContext`), field synchronisation and the thermal hooks.
"""
module Solver_Common
using TimerOutputs: @timeit
using ...Data_Manager
using ...PeriLabExceptions: @abort
using ...Helpers: check_inf_or_nan, progress_bar
using Printf: @sprintf
using ProgressBars: set_multiline_postfix, set_postfix
using ...InputDeck: SolverParams, SolverOptions, CalculationOptions
using ...ParameterSpec: model_module
import ...InputDeck: solver_name
using ...MPI_Communication: synch_responder_to_controller, synch_controller_to_responder,
                            synch_controller_bonds_to_responder_flattened
using ..Boundary_Conditions: BoundaryCondition

export TimeWindow, SolverRun, SolverContext, with_steps, synchronise_field,
       StepClock, StepFailed, prepare_solver, time_step!, time_loop, run_solver,
       create_run,
       compute_parabolic_problems_before_model_evaluation,
       compute_parabolic_problems_after_model_evaluation

"Time discretization of a solver step, computed by the solver's `init_solver`."
struct TimeWindow
    initial_time::Float64
    final_time::Float64
    dt::Float64
    nsteps::Int64
end

"""
    SolverRun{O,S}

One solver step: its parameters (`params.options isa O`), the model categories
of this step and of all steps, the stress/strain options, the time
discretization and the solver's own derived `state::S`.
"""
struct SolverRun{O<:SolverOptions,S}
    params::SolverParams
    models::Vector{String}
    all_models::Vector{String}
    calculation::CalculationOptions
    time::TimeWindow
    state::S
    function SolverRun(params::SolverParams, models::Vector{String},
                       all_models::Vector{String}, calculation::CalculationOptions,
                       time::TimeWindow, state)
        return new{typeof(params.options),typeof(state)}(params, models, all_models,
                                                         calculation, time, state)
    end
end

"""
What a solver run reads and writes besides the fields: nodes, BCs, output, and
`set_output_frequency(nsteps)`, which recomputes the output frequency when a
fallback solver takes over.
"""
struct SolverContext
    block_nodes::Dict{Int64,Vector{Int64}}
    bcs::Dict{String,BoundaryCondition}
    outputs::Dict{Int64,Dict{}}
    result_files::Vector{Dict}
    write_results::Function
    silent::Bool
    set_output_frequency::Function
end

function SolverContext(block_nodes, bcs, outputs, result_files, write_results, silent)
    return SolverContext(block_nodes, bcs, outputs, result_files, write_results, silent,
                         nsteps -> nothing)
end

"Registered name of the solver of `run`."
solver_name(run::SolverRun) = solver_name(run.params.options)

"A copy of `run` with `nsteps` steps (dry run)."
function with_steps(run::SolverRun, nsteps::Int64)
    t = run.time
    return SolverRun(run.params, run.models, run.all_models, run.calculation,
                     TimeWindow(t.initial_time, t.final_time, t.dt, nsteps), run.state)
end

"""
	synchronise_field(comm, synch_fields::Dict, overlap_map, get_field, synch_field::String, direction::String)

Synchronises field.

# Arguments
- `comm`: The MPI communicator
- `synch_fields::Dict`: A dictionary of fields
- `overlap_map`: The overlap map
- `get_field`: The function to get the field
- `synch_field::String`: The field
- `direction::String`: The direction
# Returns
- `nothing`
"""
function synchronise_field(comm,
                           synch_fields::Dict,
                           overlap_map,
                           get_field,
                           synch_field::String,
                           direction::String)
    # might not needed
    @timeit "synchronise_field" begin
        if !haskey(synch_fields, synch_field)
            @abort "Field $synch_field does not exist in synch_field dictionary"
            return nothing
        end
        if direction == "download_from_cores"
            if synch_fields[synch_field][direction]
                vector = get_field(synch_field, synch_fields[synch_field]["time"])
                @timeit "synch_responder_to_controller" return synch_responder_to_controller(comm,
                                                                                             overlap_map,
                                                                                             vector,
                                                                                             synch_fields[synch_field]["dof"])
            end
            return nothing
        end
        if direction == "upload_to_cores"
            if synch_fields[synch_field][direction]
                vector = get_field(synch_field, synch_fields[synch_field]["time"])
                if occursin("Bond", synch_field)
                    @timeit "synch_controller_bonds_to_responder_flattened" return synch_controller_bonds_to_responder_flattened(comm,
                                                                                                                                 overlap_map,
                                                                                                                                 vector,
                                                                                                                                 synch_fields[synch_field]["dof"])
                else
                    @timeit "synch_controller_to_responder" return synch_controller_to_responder(comm,
                                                                                                 overlap_map,
                                                                                                 vector,
                                                                                                 synch_fields[synch_field]["dof"])
                end
            end
            return nothing
        end
        @abort "Wrong direction key word $direction in function synchronise_field; it should be download_from_cores or upload_to_cores"
        return nothing
    end
end

"Sets the temperature of the step before the models are evaluated."
function compute_parabolic_problems_before_model_evaluation(active_nodes, run::SolverRun)
    if !("Thermal" in run.models) &&
       !("Thermal" in run.all_models)
        return
    end
    temperatureN = Data_Manager.get_field("Temperature", "N")
    temperatureNP1 = Data_Manager.get_field("Temperature", "NP1")
    deltaT = Data_Manager.get_field("Delta Temperature")
    if "Thermal" in run.models
        temperatureNP1[active_nodes] = temperatureN[active_nodes] + deltaT[active_nodes]
    else
        if "Thermal" in run.all_models
            temperatureNP1[active_nodes] = temperatureN[active_nodes]
        end
    end
end

"Computes the temperature increment from the heat flow after the models are evaluated."
function compute_parabolic_problems_after_model_evaluation(active_nodes, run::SolverRun, dt)
    if !("Thermal" in run.models) &&
       !("Thermal" in run.all_models)
        return
    end
    deltaT = Data_Manager.get_field("Delta Temperature")
    flowNP1 = Data_Manager.get_field("Heat Flow", "NP1")
    density = Data_Manager.get_field("Density")
    heat_capacity = Data_Manager.get_field("Specific Heat Capacity")
    if "Thermal" in run.models
        check_inf_or_nan(flowNP1, "Heat Flow")
        # heat capacity check. if it is zero deltaT = 0
        @views deltaT[active_nodes] = -flowNP1[active_nodes] .* dt ./
                                      (density[active_nodes] .*
                                       heat_capacity[active_nodes])
        # if fem_option && time == 0
        #     @warn "Thermal models are not supported for FEM yet."
        # end
    end
end


"Time of one step: its number, start time, time since the start of the run and increment."
struct StepClock
    idt::Int64
    time::Float64
    step_time::Float64
    dt::Float64
end

"Returned by `time_step!` when the step cannot be computed."
struct StepFailed
    reason::String
end

"""
    prepare_solver(run, ctx)

Called once per run before the first step; returns the solver's workspace
(fields, buffers, start values). Each solver adds a method.
"""
function prepare_solver end

"""
    time_step!(ws, run, ctx, clock::StepClock)

Computes one step from N to NP1. Returns `StepFailed(reason)` if the step
cannot be computed, anything else on success. Each solver adds a method.
"""
function time_step! end

"""
    run_solver(run, ctx)

Runs the solver step `run` and returns the result files. The default is
`time_loop`; a solver with its own time control may add a method.
"""
run_solver(run::SolverRun, ctx::SolverContext) = time_loop(run, ctx)

"""
    create_run(params, models, all_models, calculation, bcs, block_nodes)

Initializes the solver selected by `params.options` (its module's
`init_solver`) and returns the `SolverRun` of the step.
"""
function create_run(params::SolverParams, models::Vector{String},
                    all_models::Vector{String}, calculation::CalculationOptions,
                    bcs::Dict{String,BoundaryCondition},
                    block_nodes::Dict{Int64,Vector{Int64}})
    mod = model_module(params.options)
    # invokelatest: the solver may come from a module loaded at runtime (licensed modules)
    time, state = Base.invokelatest(mod.init_solver, params.options, params, models, bcs,
                                    block_nodes)
    return SolverRun(params, models, all_models, calculation, time, state)
end

"`p` with the given fields replaced."
function with_fields(p::T; kwargs...) where {T}
    return T((haskey(kwargs, f) ? kwargs[f] : getfield(p, f) for f in fieldnames(T))...)
end

"""
    fallback_run(run, ctx, t_N)

The run of the fallback solver of `run`, from `t_N` to the final time of `run`,
with the models, BCs and nodes of `run`.
"""
function fallback_run(run::SolverRun, ctx::SolverContext, t_N::Float64)
    # the fallback ends where the solver step ends, also if the step is given by `Additional Time`
    params = with_fields(run.params; options = run.params.fallback.options,
                         initial_time = t_N, final_time = run.time.final_time,
                         additional_time = nothing, fallback = nothing)
    return create_run(params, run.models, run.all_models, run.calculation, ctx.bcs,
                      ctx.block_nodes)
end

"""
    time_loop(run, ctx)

The time loop shared by all solvers: steps, damage limit, results, cancel,
`NP1 → N`, time and logging. A failed step is discarded; the fallback solver,
if given, recomputes it and runs to the final time, otherwise the run cancels.
"""
function time_loop(run::SolverRun, ctx::SolverContext)
    atexit(() -> Data_Manager.set_cancel(true))
    final_time = run.time.final_time
    result_files = ctx.result_files
    step_time = 0.0
    # a solver may add up external loads in a constant field; the fallback starts from the
    # loads of the step start, as if it had run the whole step
    loads = Data_Manager.has_key("External Force Densities") ?
            copy(Data_Manager.get_field("External Force Densities")) : nothing
    while true
        @info "Run $(solver_name(run)) Solver"
        result_files, failure = run_steps(run, ctx, result_files, step_time)
        failure === nothing && break
        clock, failed = failure
        @info "Step $(clock.idt) failed: $(failed.reason)"
        if run.params.fallback === nothing
            Data_Manager.set_cancel(true)
            break
        end
        @info "Switching to $(solver_name(run.params.fallback.options)) at t = $(@sprintf("%.3e", clock.time)) s"
        Data_Manager.set_current_time(clock.time)
        loads === nothing || (Data_Manager.get_field("External Force Densities") .= loads)
        step_time = clock.step_time
        run = Base.invokelatest(fallback_run, run, ctx, clock.time)
        ctx.set_output_frequency(run.time.nsteps)
    end
    Data_Manager.set_current_time(final_time)
    return result_files
end

"Maximum NP1 damage over the active nodes (0 without active nodes)."
function max_active_damage()
    damage = Data_Manager.get_damage("NP1")::Vector{Float64}
    active = Data_Manager.get_field("Active")::Vector{Bool}
    max_damage = 0.0
    nnodes = Data_Manager.get_nnodes()::Int64
    @inbounds for i in 1:nnodes
        active[i] && (max_damage = max(max_damage, damage[i]))
    end
    return max_damage
end

"""
    run_steps(run, ctx, result_files, step_time)

Steps of one run, starting at `step_time` (time since the start of the solver step). Returns the result files and `nothing`, or
`(clock, StepFailed)` of the failed step after resetting NP1.
"""
function run_steps(run::SolverRun, ctx::SolverContext, result_files, step_time::Float64)
    ws = prepare_solver(run, ctx)
    name = solver_name(run)
    rank = Data_Manager.get_rank()
    show = rank == 0 && !ctx.silent
    nsteps = run.time.nsteps
    dt = run.time.dt
    time = run.time.initial_time
    max_damage = 0.0
    damage_init = false
    damage_models = "Damage" in run.models
    iter = progress_bar(rank, nsteps, ctx.silent)
    for idt in iter
        Data_Manager.set_iteration(idt)
        clock = StepClock(idt, time, step_time, dt)
        result = @timeit name time_step!(ws, run, ctx, clock)
        if result isa StepFailed
            Data_Manager.reset_NP1()
            return result_files, (clock, result)
        end
        if damage_models
            max_damage = max_active_damage()
            if max_damage > run.params.maximum_damage
                @info "Maximum damage reached at step $idt: $max_damage"
                Data_Manager.set_cancel(true)
            end
            if !damage_init && max_damage > 0
                damage_init = true
                show && set_multiline_postfix(iter, "Damage initiated in step $idt [$time s]!")
            end
        end
        @timeit "write_results" result_files=ctx.write_results(result_files, time,
                                                               max_damage, ctx.outputs)
        if Data_Manager.get_cancel()
            show && set_multiline_postfix(iter, "Simulation canceled!")
            break
        end
        @timeit "switch_NP1_to_N" Data_Manager.switch_NP1_to_N()
        time += dt
        step_time += dt
        Data_Manager.set_current_time(time)
        if idt % ceil(nsteps / 100) == 0 || idt == 1 || idt == nsteps + 1
            @info "Step: $idt / $(nsteps+1) [$(@sprintf("%.3e", time)) s]"
        end
        show && set_postfix(iter, t = @sprintf("%.3e", time))
    end
    return result_files, nothing
end

end
