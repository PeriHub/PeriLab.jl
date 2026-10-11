# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Solver_template
using ...Data_Manager
using ...ParameterSpec: @params, register_solver
using ...InputDeck: SolverOptions, SolverParams, start_time, end_time
using ..Solver_Common: TimeWindow, SolverRun, SolverContext
import ..Solver_Common: run_solver

"""
    SolverTemplateParams

The keys of the solver's section in the input deck, e.g.
`Solver: {Final Time: 1.0, Solver Template: {Safety Factor: 0.9}}`.
"""
@params struct SolverTemplateParams <: SolverOptions
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
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
    run_solver(run, ctx)

Runs the time loop of the step `run` and returns the result files.
"""
function run_solver(run::SolverRun{SolverTemplateParams}, ctx::SolverContext)
    @info "please add your solver here"
    return ctx.result_files
end

end
