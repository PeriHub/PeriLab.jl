# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Boundary_Conditions
export init_BCs
export apply_bc_dirichlet
export apply_bc_neumann
export find_bc_free_dof
export BoundaryCondition

using ...Data_Manager
using ...InputDeck: BoundaryConditionParams, bc_node_set_names, bc_step_ids
using ...PeriLabExceptions: @abort

"""
    BoundaryCondition

A boundary condition prepared for one solver step: the resolved field
(`variable`, `time` = `"NP1"` or `"Constant"`), its `type`, the local nodes it
acts on, and its `value` — a number or an expression string from the input
deck, replaced by the compiled expression after the first evaluation.
"""
mutable struct BoundaryCondition
    variable::String
    time::String
    type::String
    initial::Bool
    coordinate::Union{Nothing,String}
    node_set::Vector{Int64}
    value::Any
end

"""
    find_bc_free_dof(bcs)

Finds all dof without a displacement boundary condition and stores them in the
Data_Manager.
"""
function find_bc_free_dof(bcs::Dict{String,BoundaryCondition})
    nnodes = Data_Manager.get_nnodes()
    dof = Data_Manager.get_dof()
    bc_free_dof = vec([(i, j) for i in 1:nnodes, j in 1:dof])
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    for bc in values(bcs)
        if bc.variable == "Displacements" && bc.type == "Dirichlet"
            act = Vector{Tuple{Int64,Int64}}([(node, dof_mapping[bc.coordinate])
                                              for node in bc.node_set])
            bc_free_dof = setdiff(bc_free_dof, act)
        end
    end
    Data_Manager.set_bc_free_dof([(t[1] + (t[2] - 1) * nnodes) for t in bc_free_dof])
end

"""
    boundary_condition(bcs_in)

Local node ids of every boundary condition (all node sets of its `Node Set`
list, in order). Aborts if a node set does not exist.
"""
function boundary_condition(bcs_in::Dict{String,BoundaryConditionParams})
    nsets = Data_Manager.get_nsets()
    node_sets = Dict{String,Vector{Int64}}()
    for (name, bc) in bcs_in
        nodes = Int64[]
        for node_set_name in bc_node_set_names(bc)
            if !haskey(nsets, node_set_name)
                @abort "Node Set '$node_set_name' is missing"
                return
            end
            append!(nodes, Data_Manager.get_local_nodes(nsets[node_set_name]))
        end
        node_sets[name] = nodes
    end
    return node_sets
end

"""
    check_valid_bcs(bcs_in, node_sets)

The boundary conditions active in the current solver step, with their field
resolved. A `z` condition in a 2D run is skipped with a warning; a missing
`Type` means Dirichlet. Aborts if the field does not exist.
"""
function check_valid_bcs(bcs_in::Dict{String,BoundaryConditionParams},
                         node_sets::Dict{String,Vector{Int64}})
    working_bcs = Dict{String,BoundaryCondition}()
    step = Data_Manager.get_step()
    dof = Data_Manager.get_dof()
    for (name, bc) in bcs_in
        steps = bc_step_ids(bc)
        if steps !== nothing && step != -1 && !(step in steps)
            continue
        end
        if bc.coordinate == "z" && dof < 3
            @warn "Boundary condition $name is not possible with $dof DOF"
            continue
        end
        type = bc.type
        if type === nothing
            type = "Dirichlet"
            @warn "Missing boundary condition type for $name. Assuming Dirichlet."
        end
        time = nothing
        for data_entry in Data_Manager.get_all_field_keys()
            if bc.variable * "NP1" == data_entry
                time = "NP1"
                break
            elseif bc.variable == data_entry
                time = "Constant"
                break
            end
        end
        if time === nothing
            @abort "Boundary condition $name is not valid: Variable $(bc.variable) not found. Please check if the physical model is activated."
            return
        end
        working_bcs[name] = BoundaryCondition(bc.variable, time, type, type == "Initial",
                                              bc.coordinate, node_sets[name], bc.value)
    end
    return working_bcs
end

"""
    init_BCs(bcs_in)

The boundary conditions of the current solver step.
"""
function init_BCs(bcs_in::Dict{String,BoundaryConditionParams})
    return check_valid_bcs(bcs_in, boundary_condition(bcs_in))
end

"""
    apply_bc_dirichlet(allowed_variables::Vector{String}, bcs::Dict{String,BoundaryCondition}, time::Float64, step_time::Float64)

Apply the boundary conditions

# Arguments
- `bcs::Dict{String,BoundaryCondition}`: The boundary conditions
- `time::Float64`: The current time
"""
function apply_bc_dirichlet(allowed_variables::Vector{String},
                            bcs::Dict{String,BoundaryCondition},
                            time::Float64,
                            step_time::Float64)
    dof = Data_Manager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = Data_Manager.get_field("Coordinates")
    for (name, bc) in bcs
        if !(bc.type in ["Initial", "Dirichlet"])
            continue
        end
        if !(bc.variable in allowed_variables)
            continue
        end
        if bc.variable == "Forces"
            field = Data_Manager.get_field("External Forces")
        elseif bc.variable == "Force Densities"
            field = Data_Manager.get_field("External Force Densities")
        else
            field = Data_Manager.get_field(bc.variable, bc.time)
        end
        if ndims(field) > 1
            if haskey(dof_mapping, bc.coordinate)
                @views field_to_apply_bc = field[bc.node_set,
                dof_mapping[bc.coordinate]]
                bc.value = eval_bc!(field_to_apply_bc,
                                       bc.value,
                                       coordinates[bc.node_set, :],
                                       time,
                                       step_time,
                                       dof,
                                       bc.initial,
                                       name)
            else
                @abort "Coordinate in boundary condition must be x,y or z."
            end
        else
            @views field_to_apply_bc = field[bc.node_set]
            bc.value = eval_bc!(field_to_apply_bc,
                                   bc.value,
                                   coordinates[bc.node_set, :],
                                   time,
                                   step_time,
                                   dof,
                                   bc.initial,
                                   name)
        end
    end
end

"""
    apply_bc_neumann(bcs::Dict{String,BoundaryCondition}, time::Float64, step_time::Float64)

Apply the boundary conditions

# Arguments
- `bcs::Dict{String,BoundaryCondition}`: The boundary conditions
- `time::Float64`: The current time
"""
function apply_bc_neumann(bcs::Dict{String,BoundaryCondition}, time::Float64,
                          step_time::Float64)
    # Currently not supported
    dof = Data_Manager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = Data_Manager.get_field("Coordinates")
    for (name, bc) in bcs
        if bc.type != "Neumann"
            continue
        end
        field = Data_Manager.get_field(bc.variable)

        if ndims(field) > 1
            if haskey(dof_mapping, bc.coordinate)
                @views field_to_apply_bc = field[bc.node_set,
                dof_mapping[bc.coordinate]]
                bc.value = eval_bc!(field_to_apply_bc,
                                       bc.value,
                                       coordinates[bc.node_set, :],
                                       time,
                                       step_time,
                                       dof,
                                       bc.initial,
                                       name,
                                       true)
            else
                @abort "Coordinate in boundary condition must be x,y or z."
                return nothing
            end
        else
            @views field_to_apply_bc = field[bc.node_set]
            bc.value = eval_bc!(field_to_apply_bc,
                                   bc.value,
                                   coordinates[bc.node_set, :],
                                   time,
                                   step_time,
                                   dof,
                                   bc.initial,
                                   name,
                                   true)
        end
    end
end

"""
    clean_up(bc::String)

Clean up the boundary condition

# Arguments
- `bc::String`: The boundary condition
# Returns
- `bc::String`: The cleaned up boundary condition
"""
function clean_up(bc::String)
    bc = replace(bc, ".*" => "*")
    bc = replace(bc, "./" => "/")
    bc = replace(bc, ".+" => "+")
    bc = replace(bc, ".-" => "-")
    bc = replace(bc, ".^" => "^")
    bc = replace(bc, ".sin" => "sin")
    bc = replace(bc, ".cos" => "cos")
    bc = replace(bc, ".tan" => "tan")
    bc = replace(bc, ".asin" => "asin")
    bc = replace(bc, ".acos" => "acos")
    bc = replace(bc, ".atan" => "atan")
    # set space before the operator to avoid integer and float problems, because the dot is connected to the number and not the operator
    bc = replace(bc, "*" => " .* ")
    bc = replace(bc, "/" => " ./ ")
    bc = replace(bc, "+" => " .+ ")
    bc = replace(bc, "-" => " .- ")
    bc = replace(bc, "^" => " .^ ")
    bc = replace(bc, "sin" => "sin.")
    bc = replace(bc, "cos" => "cos.")
    bc = replace(bc, "tan" => "tan.")
    # to guarantee the scientific number notation
    bc = replace(bc, "e .- " => "e-")
    bc = replace(bc, "e .+ " => "e+")
    bc = replace(bc, "E .- " => "e-")
    bc = replace(bc, "E .+ " => "e+")
    return bc
end

"""
    eval_bc!(field_values::Union{NodeScalarField{Float64},NodeScalarField{Int64}}, bc::Union{Float64,Float64,Int64,String}, coordinates::Matrix{Float64}, time::Float64, dof::Int64)
Working with if-statements
"if t>2 0 else 20 end"
works for scalars. If you want to evaluate a vector, please use the Julia notation as input
"ifelse.(x .> y, 10, 20)"
"""
function eval_bc!(field_values::Union{SubArray,NodeScalarField{Float64},
                                      NodeScalarField{Int64}},
                  bc::Union{Float64,Int64,String},
                  coordinates::Matrix{Float64},
                  time::Float64,
                  step_time::Float64,
                  dof::Int64,
                  initial::Bool,
                  name::String = "BC_1",
                  neumann::Bool = false)
    # reason for global
    # https://stackoverflow.com/questions/60105828/julia-local-variable-not-defined-in-expression-eval
    # the yaml input allows multiple types. But for further use this input has to be a string

    if length(coordinates) == 0
        # @warn "Ignoring boundary condition $name.\n No nodes found, check Input Deck and or Node Sets."
        return bc
    end
    bc_out = bc
    bc = string(bc)
    bc = clean_up(bc)
    if dof < 3 && occursin(r"\bz\b", bc)
        @abort "z is not valid in a 2D problem."
        return nothing
    end
    bc_value = Meta.parse(bc)

    if dof > 2
        func_args = [:x, :y, :z, :t, :st]
        dynamic_func_expr = quote
            ($(func_args...),) -> $bc_value
        end

        dynamic_bc_3D_func = Base.eval(@__MODULE__, dynamic_func_expr)

        value = Base.invokelatest(dynamic_bc_3D_func,
                                  (coordinates[:, 1], coordinates[:, 2], coordinates[:, 3],
                                   time,
                                   step_time)...)
        bc_out = dynamic_bc_3D_func
    else
        func_args = [:x, :y, :t, :st]
        dynamic_func_expr = quote
            ($(func_args...),) -> $bc_value
        end

        dynamic_2D_bc_func = Base.eval(@__MODULE__, dynamic_func_expr)

        value = Base.invokelatest(dynamic_2D_bc_func,
                                  (coordinates[:, 1], coordinates[:, 2],
                                   time,
                                   step_time)...)
        bc_out = dynamic_2D_bc_func
    end

    if isnothing(value) || (initial && time != 0.0)
        return bc_out
    end

    if value isa Number
        if neumann
            field_values .+= value
        else
            fill!(field_values, value)
        end
    else
        copyto!(field_values, value)
    end

    return bc_out
end

"""
    eval_bc!(field_values::Union{NodeScalarField{Float64},NodeScalarField{Int64}}, bc::Function, coordinates::Matrix{Float64}, time::Float64, dof::Int64)
uses the already created bc function"
"""
function eval_bc!(field_values::Union{SubArray,NodeScalarField{Float64},
                                      NodeScalarField{Int64}},
                  bc::Function,
                  coordinates::Matrix{Float64},
                  time::Float64,
                  step_time::Float64,
                  dof::Int64,
                  initial::Bool,
                  name::String = "BC_1",
                  neumann::Bool = false)
    # reason for global
    # https://stackoverflow.com/questions/60105828/julia-local-variable-not-defined-in-expression-eval
    # the yaml input allows multiple types. But for further use this input has to be a string

    if length(coordinates) == 0
        # @warn "Ignoring boundary condition $name.\n No nodes found, check Input Deck and or Node Sets."
        return bc
    end

    if dof > 2
        #func_args = [:x, :y, :z, :t, :st]
        #dynamic_func_expr = quote
        #    ($(func_args...),) -> $bc_value
        #end
        #dynamic_bc_3D_func = Base.eval(@__MODULE__, dynamic_func_expr)

        value = Base.invokelatest(bc,
                                  (coordinates[:, 1], coordinates[:, 2], coordinates[:, 3],
                                   time,
                                   step_time)...)
    else
        #func_args = [:x, :y, :t, :st]
        #dynamic_func_expr = quote
        #    ($(func_args...),) -> $bc_value
        #end
        #dynamic_2D_bc_func = Base.eval(@__MODULE__, dynamic_func_expr)
        value = Base.invokelatest(bc,
                                  (coordinates[:, 1], coordinates[:, 2],
                                   time,
                                   step_time)...)
    end

    if isnothing(value) || (initial && time != 0.0)
        return bc
    end

    if value isa Number
        if neumann
            field_values .+= value
        else
            fill!(field_values, value)
        end
    else
        copyto!(field_values, value)
    end
    return bc
end

end
