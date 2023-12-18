# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../helpers.jl")

"""
    get_FE_mesh_name(params::Dict)

Returns the name of the mesh file from the parameters

# Arguments
- `params::Dict`: The parameters
# Returns
- `String`: The name of the finite element topology file
"""
function get_FE_mesh_name(params::Dict)
    check = haskey(params["Discretization"], "Input External Topology File")
    if !check
        return nothing
    end
    return params["Discretization"]["Input External Topology File"]
end

"""
    get_mesh_name(params::Dict)

Returns the name of the mesh file from the parameters

# Arguments
- `params::Dict`: The parameters
# Returns
- `String`: The name of the mesh file
"""
function get_mesh_name(params::Dict)
    check = haskey(params["Discretization"], "Input Mesh File")
    if !check
        @error "No mesh file is defined."
        return nothing
    end
    return params["Discretization"]["Input Mesh File"]
end

"""
    get_topology_name(params::Dict)

Returns the name of the topology file from the parameters

# Arguments
- `params::Dict`: The parameters
# Returns
- `check::Bool`: Whether the topology file is defined
- `topoFile::String`: The name of the topology file
"""
function get_topology_name(params::Dict)
    check = haskey(params["Discretization"], "Input FEM Topology File")
    topoFile::String = ""
    if check
        topoFile = params["Discretization"]["Input FEM Topology File"]
    end
    return check, topoFile
end

"""
    get_bond_filters(params::Dict)

Returns the bond filters from the parameters

# Arguments
- `params::Dict`: The parameters
# Returns
- `check::Bool`: Whether the bond filters are defined
- `bfList::Dict{String,Dict{String,Any}}`: The bond filters
"""
function get_bond_filters(params::Dict)
    check = haskey(params["Discretization"], "Bond Filters")
    bfList = Dict{String,Dict{String,Any}}()
    if check
        bfList = params["Discretization"]["Bond Filters"]
    end
    return check, bfList
end

"""
    get_node_sets(params::Dict, path::String)

Returns the node sets from the parameters

# Arguments
- `params::Dict`: The parameters
- `path::String`: The path to the mesh file
# Returns
- `nsets::Dict{String,Any}`: The node sets
"""
function get_node_sets(params::Dict, path::String)
    nsets = Dict{String,Any}()
    if !haskey(params["Discretization"], "Node Sets")
        return nsets
    end
    nodesets = params["Discretization"]["Node Sets"]

    for entry in keys(nodesets)
        if nodesets[entry] isa Int64 || nodesets[entry] isa Int32
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