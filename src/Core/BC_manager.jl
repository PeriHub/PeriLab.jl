# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Boundary_Conditions
export init_BCs
export apply_bc_dirichlet
export apply_bc_neumann
export find_bc_free_dof

using ...Data_Manager
using ...Parameter_Handling: get_bc_definitions

"""
    find_bc_free_dof(bcs::Dict{String,Any})
Finds all dof without a displacement boundary condition. This tuple vector is stored in the Data_Manager.

# Arguments
- `bcs::Dict{String,Any}`: The boundary conditions
# Returns

"""
function find_bc_free_dof(bcs::Dict{Any,Any})
    nnodes = Data_Manager.get_nnodes()
    dof = Data_Manager.get_dof()
    bc_free_dof = vec([(i, j) for i in 1:nnodes, j in 1:dof])
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    for bc in values(bcs)
        if bc["Variable"] == "Displacements" && bc["Type"] == "Dirichlet"
            act = Vector{Tuple{Int64,Int64}}([(node, dof_mapping[bc["Coordinate"]])
                                              for node in bc["Node Set"]])

            bc_free_dof = setdiff(bc_free_dof, act)
        end
    end
    Data_Manager.set_bc_free_dof([(t[1] + (t[2] - 1) * nnodes) for t in bc_free_dof])
end

"""
    check_valid_bcs(bcs::Dict{String,Any})

Check if the boundary conditions are valid

# Arguments
- `bcs::Dict{String,Any}`: The boundary conditions
# Returns
- `working_bcs::Dict{String,Any}`: The valid boundary conditions
"""
function check_valid_bcs(bcs::Dict{String,Any})
    # check bc
    working_bcs = Dict()
    for bc in keys(bcs)
        if haskey(bcs[bc], "Step ID")
            if !isnothing(Data_Manager.get_step()) &&
               !(string(Data_Manager.get_step()) in split(string(bcs[bc]["Step ID"]), ","))
                continue
            end
        end

        if haskey(bcs[bc], "Coordinate")
            dof = Data_Manager.get_dof()
            if bcs[bc]["Coordinate"] == "z" && dof < 3
                @warn "Boundary condition $bc is not possible with $dof DOF"
                break
            end
        end
        valid = false

        if !haskey(bcs[bc], "Type") ||
           !(bcs[bc]["Type"] in ["Initial", "Dirichlet", "Neumann"])
            bcs[bc]["Type"] = "Dirichlet"
            @warn "Missing boundary condition type for $bc. Assuming Dirichlet."
        end
        for data_entry in Data_Manager.get_all_field_keys()
            if bcs[bc]["Variable"] * "NP1" == data_entry
                bcs[bc]["Time"] = "NP1"
                valid = true
            elseif bcs[bc]["Variable"] == data_entry
                bcs[bc]["Variable"] = data_entry
                bcs[bc]["Time"] = "Constant"
                valid = true
            end
            if valid
                bcs[bc]["Initial"] = bcs[bc]["Type"] == "Initial"
                working_bcs[bc] = bcs[bc]
                break
            end
        end
        if !valid
            @error "Boundary condition $bc is not valid: Variable $(bcs[bc]["Variable"]) not found. Please check if the physical model is activated."
            return nothing
        end
    end
    return working_bcs
end

"""
    init_BCs(params::Dict)

Initialize the boundary conditions

# Arguments
- `params::Dict`: The parameters
# Returns
- `bcs::Dict{Any,Any}`: The boundary conditions
"""
function init_BCs(params::Dict)
    bcs = boundary_condition(params)
    valid_bcs = check_valid_bcs(bcs)
    return valid_bcs
end

"""
    boundary_condition(params::Dict)

Initialize the boundary condition

# Arguments
- `params::Dict`: The parameters
# Returns
- `bcs_out::Dict{Any,Any}`: The boundary conditions
"""
function boundary_condition(params::Dict)
    bcs_in = get_bc_definitions(params)
    bcs_out = Dict{String,Any}()
    nsets = Data_Manager.get_nsets()

    for bc in keys(bcs_in)
        node_set_names = split(bcs_in[bc]["Node Set"], "+")
        node_set_names = map(r -> strip(r), node_set_names)
        bcs_out[bc] = Dict{String,Any}("Node Set" => [])
        for node_set_name in node_set_names
            if haskey(nsets, node_set_name)
                append!(bcs_out[bc]["Node Set"],
                        Data_Manager.get_local_nodes(nsets[node_set_name]))
                for entry in keys(bcs_in[bc])
                    if entry != "Node Set"
                        bcs_out[bc][entry] = bcs_in[bc][entry]
                    end
                end
            else
                @error "Node Set '$node_set_name' is missing"
                return nothing
            end
        end
    end
    return bcs_out
end

"""
    apply_bc_dirichlet(bcs::Dict, time::Float64)

Apply the boundary conditions

# Arguments
- `bcs::Dict{Any,Any}`: The boundary conditions
- `time::Float64`: The current time
"""
function apply_bc_dirichlet(allowed_variables::Vector{String},
                            bcs::Dict,
                            time::Float64,
                            step_time::Float64)
    dof = Data_Manager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = Data_Manager.get_field("Coordinates")
    for name in keys(bcs)
        bc = bcs[name]
        if !(bc["Type"] in ["Initial", "Dirichlet"])
            continue
        end
        if !(bc["Variable"] in allowed_variables)
            continue
        end
        if bc["Variable"] == "Forces"
            field = Data_Manager.get_field("External Forces")
        elseif bc["Variable"] == "Force Densities"
            field = Data_Manager.get_field("External Force Densities")
        else
            field = Data_Manager.get_field(bc["Variable"], bc["Time"])
        end
        if ndims(field) > 1
            if haskey(dof_mapping, bc["Coordinate"])
                @views field_to_apply_bc = field[bc["Node Set"],
                dof_mapping[bc["Coordinate"]]]
                bc["Value"] = eval_bc!(field_to_apply_bc,
                                       bc["Value"],
                                       coordinates[bc["Node Set"], :],
                                       time,
                                       step_time,
                                       dof,
                                       bc["Initial"],
                                       name)
            else
                @error "Coordinate in boundary condition must be x,y or z."
                return nothing
            end
        else
            @views field_to_apply_bc = field[bc["Node Set"]]
            bc["Value"] = eval_bc!(field_to_apply_bc,
                                   bc["Value"],
                                   coordinates[bc["Node Set"], :],
                                   time,
                                   step_time,
                                   dof,
                                   bc["Initial"],
                                   name)
        end
    end
end

"""
    apply_bc_neumann(bcs::Dict, time::Float64)

Apply the boundary conditions

# Arguments
- `bcs::Dict{Any,Any}`: The boundary conditions
- `time::Float64`: The current time
"""
function apply_bc_neumann(bcs::Dict, time::Float64, step_time::Float64)
    # Currently not supported
    dof = Data_Manager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = Data_Manager.get_field("Coordinates")
    for name in keys(bcs)
        bc = bcs[name]
        if bc["Type"] != "Neumann"
            continue
        end
        field = Data_Manager.get_field(bc["Variable"])

        if ndims(field) > 1
            if haskey(dof_mapping, bc["Coordinate"])
                @views field_to_apply_bc = field[bc["Node Set"],
                dof_mapping[bc["Coordinate"]]]
                bc["Value"] = eval_bc!(field_to_apply_bc,
                                       bc["Value"],
                                       coordinates[bc["Node Set"], :],
                                       time,
                                       step_time,
                                       dof,
                                       bc["Initial"],
                                       name,
                                       true)
            else
                @error "Coordinate in boundary condition must be x,y or z."
                return nothing
            end
        else
            @views field_to_apply_bc = field[bc["Node Set"]]
            bc["Value"] = eval_bc!(field_to_apply_bc,
                                   bc["Value"],
                                   coordinates[bc["Node Set"], :],
                                   time,
                                   step_time,
                                   dof,
                                   bc["Initial"],
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
    if dof < 2 && "z" in bc
        @error "z is not valid in a 2D problem."
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
