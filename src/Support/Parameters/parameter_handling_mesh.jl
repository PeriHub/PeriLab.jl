# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_mesh_name(params)
    check = check_element(params["Discretization"], "Input Mesh File")
    if !check
        @error "No mesh file is defined."
        return
    end
    return params["Discretization"]["Input Mesh File"]
end

function get_topology_name(params)
    check = check_element(params["Discretization"], "Input FEM Topology File")
    topoFile::String = ""
    if check
        topoFile = params["Discretization"]["Input FEM Topology File"]
    end
    return check, topoFile
end

function get_bond_filters(params)
    check = check_element(params["Discretization"], "Bond Filters")
    bfList = Any
    if check
        bfList = params["Discretization"]["Bond Filters"]
    end
    return check, bfList
end

function get_node_sets(params)
    nsets = Dict{String,Any}()
    if check_element(params["Discretization"], "Node Sets") == false
        return []
    end
    nodesets = params["Discretization"]["Node Sets"]
    for entry in keys(nodesets)
        if (typeof(nodesets[entry]) == Int64) | (typeof(nodesets[entry]) == Int32)
            nsets[entry] = [nodesets[entry]]
        elseif occursin(".txt", nodesets[entry])
            nodes = CSV.read(nodesets[entry], DataFrame; delim=" ", header=false)
            if size(nodes) == (0, 0)
                @error "Node set file is empty " * nodesets[entry]
            end
            nsets[entry] = nodes.Column1
        else
            nodes = split(nodesets[entry])
            nsets[entry] = parse.(Int, nodes)
        end
    end
    return nsets
end