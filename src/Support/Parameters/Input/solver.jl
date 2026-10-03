# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Every solver option section carries "Safety Factor", "Fixed dt" and
# "Numerical Damping": the solver getters read them from the active section.

@params struct VerletParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
end

@params struct StaticParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    solution_tolerance::Float64 = opt("Solution tolerance"; default = 1e-7, min = 0)
    residual_tolerance::Float64 = opt("Residual tolerance"; default = 1e-7, min = 0)
    maximum_number_of_iterations::Int64 = opt("Maximum number of iterations"; default = 100,
                                              min = 1)
    show_solver_iteration::Bool = opt("Show solver iteration"; default = false)
    residual_scaling::Float64 = opt("Residual scaling"; default = 1e6)
    m::Int64 = opt("m"; default = 15, min = 1)
    linear_start_value::Union{Nothing,String} = opt("Linear Start Value"; default = nothing)
    nlsolve::Union{Nothing,Bool} = opt("NLSolve"; default = nothing)
    solver_type::Union{Nothing,String} = opt("Solver Type"; default = nothing)
end

@params struct LinearStaticMatrixParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    matrix_update::Bool = opt("Matrix Update"; default = false)
end

@params struct ModelReductionParams
    type::String = req("Type")
    number_of_modes::Int64 = opt("Number of Modes"; default = 1, min = 1)
    material_point_region::Bool = opt("Material Point Region"; default = true)
    reduction_blocks::Union{Nothing,String} = opt("Reduction Blocks"; default = nothing)
end

@params struct VerletMatrixParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    model_reduction::Union{Nothing,ModelReductionParams} = opt("Model Reduction";
                                                               default = nothing)
end

@params struct NewmarkParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    matrix_update::Bool = opt("Matrix Update"; default = false)
    newmark_delta::Float64 = opt("Newmark Delta"; default = 0.5)
    newmark_alpha::Union{Nothing,Float64} = opt("Newmark Alpha"; default = nothing)
end

"Used for `Solver` and for every step of `Multistep Solver`."
@params struct SolverParams
    initial_time::Union{Nothing,Float64} = opt("Initial Time"; default = nothing,
                                               quantity = :time)
    final_time::Union{Nothing,Float64} = opt("Final Time"; default = nothing, quantity = :time)
    additional_time::Union{Nothing,Float64} = opt("Additional Time"; default = nothing,
                                                  quantity = :time)
    number_of_steps::Int64 = opt("Number of Steps"; default = 1, min = 1)
    maximum_damage::Float64 = opt("Maximum Damage"; default = Inf)
    step_id::Union{Nothing,Int64} = opt("Step ID"; default = nothing)
    additive_models::Bool = opt("Additive Models"; default = false)
    degradation_models::Bool = opt("Degradation Models"; default = false)
    damage_models::Bool = opt("Damage Models"; default = false)
    material_models::Bool = opt("Material Models"; default = true)
    thermal_models::Bool = opt("Thermal Models"; default = false)
    pre_calculation_models::Bool = opt("Pre Calculation Models"; default = true)
    calculate_cauchy::Bool = opt("Calculate Cauchy"; default = false)
    calculate_von_mises_stress::Bool = opt("Calculate von Mises stress"; default = false)
    calculate_strain::Bool = opt("Calculate Strain"; default = false)
    verlet::Union{Nothing,VerletParams} = opt("Verlet"; default = nothing)
    static::Union{Nothing,StaticParams} = opt("Static"; default = nothing)
    linear_static_matrix_based::Union{Nothing,LinearStaticMatrixParams} = opt("Linear Static Matrix Based";
                                                                              default = nothing)
    verlet_matrix_based::Union{Nothing,VerletMatrixParams} = opt("Verlet Matrix Based";
                                                                 default = nothing)
    newmark::Union{Nothing,NewmarkParams} = opt("Newmark"; default = nothing)
end

const SOLVER_NAMES = ("Verlet", "Static", "Linear Static Matrix Based", "Verlet Matrix Based",
                      "Newmark")

function check!(p::SolverParams, path::String, ctx::ParseContext)
    sections = (p.verlet, p.static, p.linear_static_matrix_based, p.verlet_matrix_based,
                p.newmark)
    given = [name for (name, section) in zip(SOLVER_NAMES, sections) if section !== nothing]
    if isempty(given)
        add_error!(ctx, path,
                   "one solver is required: Verlet, Static, Linear Static Matrix Based, Verlet Matrix Based or Newmark")
    elseif length(given) > 1
        add_error!(ctx, path, "only one solver may be given, found: $(join(given, ", "))")
    end
    if p.final_time === nothing && p.additional_time === nothing
        add_error!(ctx, path, "\"Final Time\" or \"Additional Time\" is required")
    end
    return nothing
end
