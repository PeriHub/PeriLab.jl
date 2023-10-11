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
        if bcs[bc]["Coordinate"] == "z" && datamanager.get_dof() < 3
            break
        end
        for dataentry in datamanager.get_all_field_keys()
            if (occursin(bcs[bc]["Type"], dataentry) && occursin("NP1", dataentry)) || bcs[bc]["Type"] == dataentry
                working_bcs[bc] = bcs[bc]
                bcs[bc]["Type"] = dataentry
                break
            end
        end
    end
    return working_bcs
end

function init_BCs(params, datamanager)
    bcs = boundary_condition(params, datamanager)
    bcs = check_valid_bcs(bcs, datamanager)
    return bcs
end

function boundary_condition(params, datamanager)
    bcs_in = get_bc_definitions(params)
    bcs_out = Dict{String,Any}()
    nsets = datamanager.get_nsets()

    for bc in keys(bcs_in)
        node_set_name = bcs_in[bc]["Node Set"]
        bcs_out[bc] = Dict{String,Any}("Node Set" => datamanager.get_local_nodes(nsets[node_set_name]))
        for entry in keys(bcs_in[bc])
            if entry != "Node Set"
                bcs_out[bc][entry] = bcs_in[bc][entry]
            end
        end
    end
    return bcs_out
end

function apply_bc(bcs, datamanager, time)
    dof = datamanager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = datamanager.get_field("Coordinates")
    for name in keys(bcs)
        bc = bcs[name]
        field_to_apply_bc = datamanager.get_field(bc["Type"])
        if length(field_to_apply_bc) == 0
            field_to_apply_bc = datamanager.get_field(bc["Type"])
        end
        field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]] = eval_bc(bc["Value"], coordinates[bc["Node Set"], :], time, dof)
    end
end

function clean_up(bc::String)
    bc = replace(bc, "*" => " .* ")
    bc = replace(bc, "/" => " ./ ")
    bc = replace(bc, "+" => " .+ ")
    bc = replace(bc, "-" => " .- ")
    return bc
end

function eval_bc(bc::Union{Float32,Float64,Int64,String}, coordinates, time, dof)
    # reason for global
    # https://stackoverflow.com/questions/60105828/julia-local-variable-not-defined-in-expression-eval
    # the yaml input allows multiple types. But for further use this input has to be a string
    bc = string(bc)
    bc = clean_up(bc)
    bc_value = Meta.parse(bc)
    """
    Working with if-statements
      "if t>2 0 else 20 end"
      works for scalars. If you want to evaluate a vector, please use the Julia notation as input
      "ifelse.(x .> y, 10, 20)"
    """
    global x = coordinates[:, 1]
    global y = coordinates[:, 2]
    global t = time

    if dof > 2
        global z = coordinates[:, 3]
    else
        global z = zeros(typeof(x[1]), length(x))
    end
    return zeros(Float32, length(x)) .+ eval(bc_value)
end

end