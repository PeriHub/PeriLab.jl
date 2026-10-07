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

A validated input deck. `materials` holds the typed material models
(`ParameterSpec.WithBase`), keyed by name; `damages`, `additives` and `degradations`
hold the typed damage, additive and degradation models;
`models` stays the raw `Models` dict for the categories not yet migrated;
`globals` is the unvalidated `Globals` escape hatch.
"""
struct PeriLabInput
    sections::PeriLabSections
    contact::Union{Nothing,ContactInput}
    materials::Dict{String,Any}
    damages::Dict{String,Any}
    additives::Dict{String,Any}
    degradations::Dict{String,Any}
    models::Dict{String,Any}
    globals::Dict{String,Any}
end

const _SPECIAL_KEYS = ("Models", "Contact", "Globals")

"Typed models of `section` (e.g. \"Material Models\"), keyed by name; errors go to `ctx`."
function parse_models(models::AbstractDict, section::String, category::Symbol,
                      name_key::String, ctx::ParseContext)
    parsed = Dict{String,Any}()
    raw = get(models, section, nothing)
    raw === nothing && return parsed
    path = join_path("Models", section)
    if !(raw isa AbstractDict)
        add_error!(ctx, path,
                   "expected a section of `key: value` entries, got $(ParameterSpec._describe(raw))")
        return parsed
    end
    for (name, entry) in raw
        entry_path = join_path(path, string(name))
        if !(entry isa AbstractDict)
            add_error!(ctx, entry_path,
                       "expected a section of `key: value` entries, got $(ParameterSpec._describe(entry))")
            continue
        end
        model = ParameterSpec.parse_model(category,
                                          Dict{String,Any}(string(k) => v for (k, v) in entry),
                                          entry_path, ctx; name_key = name_key)
        model === nothing || (parsed[string(name)] = model)
    end
    return parsed
end

parse_materials(models::AbstractDict, ctx::ParseContext) = parse_models(models,
                                                                        "Material Models",
                                                                        :material,
                                                                        "Material Model", ctx)

"Typed models of `section` that cannot be combined with `+`."
function parse_single_models(models::AbstractDict, section::String, category::Symbol,
                             name_key::String, ctx::ParseContext)
    parsed = parse_models(models, section, category, name_key, ctx)
    for (name, m) in parsed
        part = m isa ParameterSpec.WithBase ? m.model : m
        part isa ParameterSpec.Composite || continue
        add_error!(ctx, join_path(join_path(join_path("Models", section), name), name_key),
                   "$category models cannot be combined with +")
    end
    return parsed
end

parse_damages(models::AbstractDict, ctx::ParseContext) = parse_single_models(models,
                                                                             "Damage Models",
                                                                             :damage,
                                                                             "Damage Model",
                                                                             ctx)

"""
    read_input(deck, directory = ""; strict = true) -> (input, ctx)

Reads the content of the `PeriLab:` key of an input deck. `directory` is the
deck's directory (for relative data files). All problems are collected in
`ctx`; `input` is `nothing` if any of them is an error. Does not abort —
call `ParameterSpec.report!(ctx)` for that.
"""
function read_input(deck::AbstractDict, directory::AbstractString = ""; strict::Bool = true)
    ctx = ParseContext(directory = directory, strict = strict)
    # the special keys are read below; passing them as known keys makes
    # misspellings of them ("contact", "Contacts") errors with suggestions
    sections = parse_section(PeriLabSections, Dict{String,Any}(string(k) => v for (k, v) in deck),
                             "", ctx; known_extra = _SPECIAL_KEYS)
    contact = haskey(deck, "Contact") ? parse_contact(deck["Contact"], "Contact", ctx) :
              nothing
    models = get(deck, "Models", nothing)
    if models === nothing
        add_error!(ctx, "Models", "missing (required)")
    elseif !(models isa AbstractDict)
        add_error!(ctx, "Models",
                   "expected a section of `key: value` entries, got $(ParameterSpec._describe(models))")
    end
    materials = models isa AbstractDict ? parse_materials(models, ctx) : Dict{String,Any}()
    damages = models isa AbstractDict ? parse_damages(models, ctx) : Dict{String,Any}()
    additives = models isa AbstractDict ?
                parse_single_models(models, "Additive Models", :additive, "Additive Model",
                                    ctx) : Dict{String,Any}()
    degradations = models isa AbstractDict ?
                   parse_single_models(models, "Degradation Models", :degradation,
                                       "Degradation Model", ctx) : Dict{String,Any}()
    globals = get(deck, "Globals", Dict{String,Any}())
    if ParameterSpec.has_errors(ctx) || sections === nothing
        return nothing, ctx
    end
    input = PeriLabInput(sections, contact, materials, damages, additives, degradations,
                         Dict{String,Any}(string(k) => v for (k, v) in models),
                         globals isa AbstractDict ?
                         Dict{String,Any}(string(k) => v for (k, v) in globals) :
                         Dict{String,Any}())
    return input, ctx
end

"Step IDs to run, sorted; `[-1]` for a single `Solver`."
function solver_steps(input::PeriLabInput)
    isempty(input.sections.multistep_solver) && return [-1]
    return sort!([step.step_id for step in values(input.sections.multistep_solver)])
end

"Solver parameters of step `step_id` (`-1`: the single `Solver`)."
function solver_step(input::PeriLabInput, step_id::Int64)
    step_id == -1 && return input.sections.solver
    for step in values(input.sections.multistep_solver)
        step.step_id == step_id && return step
    end
    @abort "Step ID $step_id not found"
end
