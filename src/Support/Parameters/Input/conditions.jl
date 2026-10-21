# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct BoundaryConditionParams
    variable::String = req("Variable")
    node_set::String = req("Node Set"; description = "Node set name; several joined with +")
    value::Union{Float64,String} = req("Value"; description = "Number or expression, e.g. \"10*t\"")
    type::Union{Nothing,String} = opt("Type"; default = nothing,
                                      allowed = ["Initial", "Dirichlet", "Neumann"])
    coordinate::Union{Nothing,String} = opt("Coordinate"; default = nothing)
    step_id::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
end

@params struct ComputeClassParams
    compute_class::String = req("Compute Class")
    variable::String = req("Variable")
    calculation_type::Union{Nothing,String} = opt("Calculation Type"; default = nothing)
    block::Union{Nothing,String} = opt("Block"; default = nothing)
    node_set::Union{Nothing,String} = opt("Node Set"; default = nothing)
    equation::Union{Nothing,String} = opt("Equation"; default = nothing)
    x::Union{Nothing,Float64} = opt("X"; default = nothing)
    y::Union{Nothing,Float64} = opt("Y"; default = nothing)
    z::Union{Nothing,Float64} = opt("Z"; default = nothing)
end

function check!(bc::BoundaryConditionParams, path::String, ctx::ParseContext)
    if bc.step_id isa String &&
       any(part -> tryparse(Int64, strip(part)) === nothing, split(bc.step_id, ","))
        add_error!(ctx, join_path(path, "Step ID"),
                   "expected an integer or a comma-separated list of integers, got \"$(bc.step_id)\"")
    end
    return nothing
end

"Names of the node sets of a boundary condition (`Node Set` joined with `+`)."
bc_node_set_names(bc::BoundaryConditionParams) = String.(strip.(split(bc.node_set, "+")))

"Solver steps a boundary condition applies to, or `nothing` (all steps)."
function bc_step_ids(bc::BoundaryConditionParams)
    bc.step_id === nothing && return nothing
    bc.step_id isa Int64 && return [bc.step_id]
    return parse.(Int64, strip.(split(bc.step_id, ",")))
end

"Names of the compute classes, sorted."
compute_names(computes::Dict{String,ComputeClassParams}) = sort!(collect(keys(computes)))

"""
    active_computes(computes, field_keys)

The compute classes whose `Variable` is an existing field (directly or as
`<Variable>NP1`); the others are skipped with a warning.
"""
function active_computes(computes::Dict{String,ComputeClassParams},
                         field_keys::Vector{String})
    active = Dict{String,ComputeClassParams}()
    for (name, compute) in computes
        if compute.variable in field_keys || compute.variable * "NP1" in field_keys
            active[name] = compute
        else
            @warn '"' * compute.variable * '"' * " is not defined as variable"
        end
    end
    return active
end
