# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params "Prescribes a field, e.g. displacements, forces or temperature, on node sets" struct BoundaryConditionParams
    "Field to prescribe, e.g. \"Displacements\", \"Forces\", \"Force Densities\" or \"Temperature\""
    variable::String = req("Variable")
    "Node set name from Discretization; several joined with +"
    node_set::String = req("Node Set")
    "Number or expression in t (time), st (time within the solver step) and x, y, z (coordinates), e.g. \"10*t\""
    value::Union{Float64,String} = req("Value"; evaluated = true)
    "Initial: set once at the start; Dirichlet: set in every step; Neumann: added to the field in every step. Dirichlet if not given"
    type::Union{Nothing,String} = opt("Type"; default = nothing,
                                      allowed = ["Initial", "Dirichlet", "Neumann"])
    "Component of a vector field the condition acts on: x, y or z"
    coordinate::Union{Nothing,String} = opt("Coordinate"; default = nothing)
    "Solver steps of a Multistep Solver the condition is active in: one id or a comma-separated list; all steps if not given"
    step_id::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
end

@params "A global value computed from a field over a block, a node set or the node nearest to a point, for the outputs" struct ComputeClassParams
    "Where the value is taken: \"Block_Data\", \"Node_Set_Data\" or \"Nearest_Point_Data\""
    compute_class::String = req("Compute Class")
    "Field the value is computed from, e.g. \"Forces\" or \"Displacements\""
    variable::String = req("Variable")
    "Reduction over the nodes: \"Sum\", \"Maximum\", \"Minimum\" or \"Average\"; without it a node set exports every node separately"
    calculation_type::Union{Nothing,String} = opt("Calculation Type"; default = nothing)
    "Block name (Block_Data)"
    block::Union{Nothing,String} = opt("Block"; default = nothing)
    "Node set name (Node_Set_Data)"
    node_set::Union{Nothing,String} = opt("Node Set"; default = nothing)
    "Expression in x applied to the computed value x, e.g. \"-x\" or \"x/2\""
    equation::Union{Nothing,String} = opt("Equation"; default = nothing, evaluated = true)
    "x coordinate of the point (Nearest_Point_Data)"
    x::Union{Nothing,Float64} = opt("X"; default = nothing)
    "y coordinate of the point (Nearest_Point_Data)"
    y::Union{Nothing,Float64} = opt("Y"; default = nothing)
    "z coordinate of the point (Nearest_Point_Data, 3D)"
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
