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
using ...Helpers: check_inf_or_nan
using ...InputDeck: SolverParams, SolverOptions, CalculationOptions
import ...InputDeck: solver_name
using ...MPI_Communication: synch_responder_to_controller, synch_controller_to_responder,
                            synch_controller_bonds_to_responder_flattened
using ..Boundary_Conditions: BoundaryCondition

export TimeWindow, SolverRun, SolverContext, with_steps, synchronise_field,
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

"What a solver run reads and writes besides the fields: nodes, BCs and output."
struct SolverContext
    block_nodes::Dict{Int64,Vector{Int64}}
    bcs::Dict{String,BoundaryCondition}
    outputs::Dict{Int64,Dict{}}
    result_files::Vector{Dict}
    write_results::Function
    silent::Bool
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

end
