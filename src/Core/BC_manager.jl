# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Boundary_conditions
export init_BCs
export apply_bc
include("../Support/Parameters/parameter_handling.jl")
function check_valid_bcs(bcs, datamanager)
    # check bc
    working_bcs = Dict()
    for bc in keys(bcs)
        if haskey(bcs[bc], "Coordinate")
            dof = datamanager.get_dof()
            if bcs[bc]["Coordinate"] == "z" && dof < 3
                @warn "Boundary condition $bc is not possible with $dof DOF"
                break
            end
        end
        valid = false
        for dataentry in datamanager.get_all_field_keys()
            initial = occursin("Initial", bcs[bc]["Type"])
            bc_type = replace(bcs[bc]["Type"], "Initial " => "")
            if (occursin(bc_type, dataentry) && occursin("NP1", dataentry)) || bc_type == dataentry
                working_bcs[bc] = bcs[bc]
                bcs[bc]["Type"] = dataentry
                bcs[bc]["Initial"] = initial
                valid = true
                break
            end
        end
        if !valid
            @warn "Boundary condition $bc is not valid"
        end
    end
    return working_bcs
end

function init_BCs(params::Dict, datamanager)
    bcs = boundary_condition(params, datamanager)
    bcs = check_valid_bcs(bcs, datamanager)
    return bcs
end

function boundary_condition(params::Dict, datamanager)
    bcs_in = get_bc_definitions(params)
    bcs_out = Dict{String,Any}()
    nsets = datamanager.get_nsets()

    for bc in keys(bcs_in)
        node_set_name = bcs_in[bc]["Node Set"]
        if haskey(nsets, node_set_name)
            bcs_out[bc] = Dict{String,Any}("Node Set" => datamanager.get_local_nodes(nsets[node_set_name]))
            for entry in keys(bcs_in[bc])
                if entry != "Node Set"
                    bcs_out[bc][entry] = bcs_in[bc][entry]
                end
            end
        else
            @error "Node Set '$node_set_name' is missing"
        end
    end
    return bcs_out
end

function apply_bc(bcs::Dict, datamanager::Module, time::Float64)
    dof = datamanager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = datamanager.get_field("Coordinates")
    for name in keys(bcs)
        bc = bcs[name]
        field_to_apply_bc = datamanager.get_field(bc["Type"])
        if length(field_to_apply_bc) == 0
            field_to_apply_bc = datamanager.get_field(bc["Type"])
        end

        if ndims(field_to_apply_bc) > 1
            if haskey(dof_mapping, bc["Coordinate"])
                field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]] = eval_bc(field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]], bc["Value"], coordinates[bc["Node Set"], :], time, dof, bc["Initial"], name)
            else
                @error "Coordinate in boundary condition must be x,y or z."
                return nothing
            end
        else
            field_to_apply_bc[bc["Node Set"]] = eval_bc(field_to_apply_bc[bc["Node Set"]], bc["Value"], coordinates[bc["Node Set"], :], time, dof, bc["Initial"], name)
        end

    end
    return datamanager
end

function clean_up(bc::String)
    # set space before the operator to avoid integer and float problems, because the dot is connected to the number and not the operator
    bc = replace(bc, "*" => " .* ")
    bc = replace(bc, "/" => " ./ ")
    bc = replace(bc, "+" => " .+ ")
    bc = replace(bc, "-" => " .- ")
    # to guarantee the scientific number notation
    bc = replace(bc, "e .- " => "e-")
    bc = replace(bc, "e .+ " => "e+")
    bc = replace(bc, "E .- " => "e-")
    bc = replace(bc, "E .+ " => "e+")
    return bc
end
"""
eval_bc(field_values::Union{SubArray,Vector{Float64},Vector{Int64}}, bc::Union{Float64,Float64,Int64,String}, coordinates::Matrix{Float64}, time::Float64, dof::Int64)
Working with if-statements
  "if t>2 0 else 20 end"
  works for scalars. If you want to evaluate a vector, please use the Julia notation as input
  "ifelse.(x .> y, 10, 20)"
"""
function eval_bc(field_values::Union{SubArray,Vector{Float64},Vector{Int64}}, bc::Union{Float64,Int64,String}, coordinates::Union{Matrix{Float64},Matrix{Int64}}, time::Float64, dof::Int64, initial::Bool, name::String="BC_1")
    # reason for global
    # https://stackoverflow.com/questions/60105828/julia-local-variable-not-defined-in-expression-eval
    # the yaml input allows multiple types. But for further use this input has to be a string
    println(bc)
    bc = string(bc)
    println(bc)
    bc = clean_up(bc)
    println(bc)
    bc_value = Meta.parse(bc)

    if length(coordinates) == 0
        @warn "Ignoring boundary condition $name.\n No nodes found, check Input Deck and or Node Sets."
        return field_values
    end

    global x = coordinates[:, 1]
    global y = coordinates[:, 2]
    global t = time

    if dof > 2
        global z = coordinates[:, 3]
    else
        global z = zeros(typeof(x[1]), length(x))
    end
    println(bc, bc_value)
    value = eval(bc_value)

    if isnothing(value) || (initial && t != 0.0)
        return field_values
    end

    return zeros(Float64, length(x)) .+ value
end

end