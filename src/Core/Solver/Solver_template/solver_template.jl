# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Solver_template
using ...Data_Manager
using ...ParameterSpec: @params, register_solver
using ...InputDeck: SolverOptions, SolverParams, start_time, end_time
using ..Solver_Common: TimeWindow, SolverRun, SolverContext, StepClock, StepFailed
import ..Solver_Common: prepare_solver, time_step!

"""
    SolverTemplateParams

The keys of the solver's section in the input deck, e.g.
`Solver: {Final Time: 1.0, Solver Template: {Safety Factor: 0.9}}`.
"""
@params "Template of a solver; copy it to add a new solver" struct SolverTemplateParams <: SolverOptions
    "Factor on the stable time step; values below 1 add a safety margin"
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    "Time step to use instead of the stable time step; -1: compute it from the stable time step and the safety factor"
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    "Damps the velocity: the velocity of the previous step is scaled by one minus this value"
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
end

# Uncomment to make the solver available as `Solver Template:` in the input deck.
# __init__() = register_solver("Solver Template", SolverTemplateParams)

"""
    init_solver(options, params, models, bcs, block_nodes)

Creates the solver's fields and returns the time discretization and the
solver's derived state (`nothing` if it has none).
"""
function init_solver(options::SolverTemplateParams, params::SolverParams,
                     models::Vector{String}, bcs, block_nodes)
    initial_time = start_time(params, Data_Manager.get_current_time())
    final_time = end_time(params, Data_Manager.get_current_time())
    nsteps = something(params.number_of_steps, 1)
    return TimeWindow(initial_time, final_time, (final_time - initial_time) / nsteps, nsteps),
           nothing
end

"""
    SolverTemplateWorkspace

What the solver keeps between steps: fields fetched once, buffers, start values.
"""
struct SolverTemplateWorkspace
    volume::Vector{Float64}
end

"""
    prepare_solver(run, ctx)

Called once before the first step, after the models are initialized.
"""
function prepare_solver(run::SolverRun{SolverTemplateParams}, ctx::SolverContext)
    return SolverTemplateWorkspace(Data_Manager.get_field("Volume"))
end

"""
    time_step!(ws, run, ctx, clock)

Computes one step from N to NP1 at `clock.time` with increment `clock.dt`.
Return `StepFailed(reason)` if the step cannot be computed; the shared loop then
discards it and hands it to the `Fallback` solver of the deck, if any. Writing
results, the damage limit, cancelling, `NP1 → N` and the time are done by the
shared loop.
"""
function time_step!(ws::SolverTemplateWorkspace, run::SolverRun{SolverTemplateParams},
                    ctx::SolverContext, clock::StepClock)
    @info "please add your solver step here"
    return nothing
end

end
