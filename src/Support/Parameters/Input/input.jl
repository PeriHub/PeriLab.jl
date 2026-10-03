# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct PeriLabSections
    discretization::DiscretizationParams = req("Discretization")
    blocks::Dict{String,BlockParams} = req("Blocks")
    fem::Union{Nothing,FEMParams} = opt("FEM"; default = nothing)
    surface_correction::Union{Nothing,SurfaceCorrectionParams} = opt("Surface Correction";
                                                                     default = nothing)
    solver::Union{Nothing,SolverParams} = opt("Solver"; default = nothing)
    multistep_solver::Dict{String,SolverParams} = opt("Multistep Solver";
                                                      default = Dict{String,SolverParams}())
    outputs::Dict{String,OutputParams} = opt("Outputs"; default = Dict{String,OutputParams}())
    boundary_conditions::Dict{String,BoundaryConditionParams} = opt("Boundary Conditions";
                                                                    default = Dict{String,
                                                                                   BoundaryConditionParams}())
    compute_class_parameters::Dict{String,ComputeClassParams} = opt("Compute Class Parameters";
                                                                    default = Dict{String,
                                                                                   ComputeClassParams}())
    strict_validation::Bool = opt("Strict Validation"; default = true,
                                  description = "false reports unknown keys as warnings")
end

function check!(p::PeriLabSections, path::String, ctx::ParseContext)
    isempty(p.blocks) && add_error!(ctx, join_path(path, "Blocks"),
                                    "at least one block is required")
    if p.solver === nothing && isempty(p.multistep_solver)
        add_error!(ctx, path, "\"Solver\" or \"Multistep Solver\" is required")
    end
    if p.solver !== nothing
        solver_path = join_path(path, "Solver")
        for (key, value) in (("Initial Time", p.solver.initial_time),
                             ("Final Time", p.solver.final_time))
            value === nothing &&
                add_error!(ctx, join_path(solver_path, key), "missing (required for \"Solver\")")
        end
    end
    for (name, step) in p.multistep_solver
        if step.step_id === nothing
            add_error!(ctx,
                       join_path(join_path(join_path(path, "Multistep Solver"), name), "Step ID"),
                       "missing (required in a multistep solver step)")
        end
    end
    return nothing
end

"""
    PeriLabInput

A validated input deck. `models` stays the raw `Models` dict until the model
modules declare their parameters (phase 3); `globals` is the unvalidated
`Globals` escape hatch.
"""
struct PeriLabInput
    sections::PeriLabSections
    contact::Union{Nothing,ContactInput}
    models::Dict{String,Any}
    globals::Dict{String,Any}
end

const _SPECIAL_KEYS = ("Models", "Contact", "Globals")

"""
    read_input(deck, directory = ""; strict = true) -> (input, ctx)

Reads the content of the `PeriLab:` key of an input deck. `directory` is the
deck's directory (for relative data files). All problems are collected in
`ctx`; `input` is `nothing` if any of them is an error. Does not abort —
call `ParameterSpec.report!(ctx)` for that.
"""
function read_input(deck::AbstractDict, directory::AbstractString = ""; strict::Bool = true)
    ctx = ParseContext(directory = directory, strict = strict)
    plain = Dict{String,Any}(string(k) => v for (k, v) in deck
                             if !(string(k) in _SPECIAL_KEYS))
    sections = parse_section(PeriLabSections, plain, "", ctx)
    contact = haskey(deck, "Contact") ? parse_contact(deck["Contact"], "Contact", ctx) :
              nothing
    models = get(deck, "Models", nothing)
    if models === nothing
        add_error!(ctx, "Models", "missing (required)")
    elseif !(models isa AbstractDict)
        add_error!(ctx, "Models",
                   "expected a section of `key: value` entries, got $(ParameterSpec._describe(models))")
    end
    globals = get(deck, "Globals", Dict{String,Any}())
    if ParameterSpec.has_errors(ctx) || sections === nothing
        return nothing, ctx
    end
    input = PeriLabInput(sections, contact, Dict{String,Any}(string(k) => v for (k, v) in models),
                         globals isa AbstractDict ?
                         Dict{String,Any}(string(k) => v for (k, v) in globals) :
                         Dict{String,Any}())
    return input, ctx
end
