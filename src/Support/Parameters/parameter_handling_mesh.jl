# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../helpers.jl")

function get_mesh_name(params::Dict)
    check = check_element(params["Discretization"], "Input Mesh File")
    if !check
        @error "No mesh file is defined."
        return nothing
    end
    return params["Discretization"]["Input Mesh File"]
end

function get_topology_name(params::Dict)
    check = check_element(params["Discretization"], "Input FEM Topology File")
    topoFile::String = ""
    if check
        topoFile = params["Discretization"]["Input FEM Topology File"]
    end
    return check, topoFile
end

function get_bond_filters(params::Dict)
    check = check_element(params["Discretization"], "Bond Filters")
    bfList = Dict{String,Dict{String,Any}}()
    if check
        bfList = params["Discretization"]["Bond Filters"]
    end
    return check, bfList
end

function get_node_sets(params::Dict, path::String)
    nsets = Dict{String,Any}()
    if !check_element(params["Discretization"], "Node Sets")
        return nsets
    end
    nodesets = params["Discretization"]["Node Sets"]

    for entry in keys(nodesets)
        if (typeof(nodesets[entry]) == Int64) | (typeof(nodesets[entry]) == Int32)
            nsets[entry] = [nodesets[entry]]
        elseif occursin(".txt", nodesets[entry])

            if isnothing(get_header(joinpath(path, nodesets[entry])))
                @warn "Node set file " * nodesets[entry] * " is not correctly specified. Please check the examples. The node set is excluded."
                continue
            end
            header_line, header = get_header(joinpath(path, nodesets[entry]))
            nodes = CSV.read(joinpath(path, nodesets[entry]), DataFrame; delim=" ", header=false, skipto=header_line + 1)
            if size(nodes) == (0, 0)
                @error "Node set file is empty " * joinpath(path, nodesets[entry]) * ". The node set is excluded."
                continue
            end
            nsets[entry] = nodes.Column1
        else
            nodes = split(nodesets[entry])
            nsets[entry] = parse.(Int, nodes)
        end
    end
    return nsets
end