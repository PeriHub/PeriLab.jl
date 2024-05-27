# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
export get_mesh_name
export get_external_topology_name
export get_bond_filters
export get_node_sets
export get_header
using AbaqusReader
using Exodus
"""
    get_external_topology_name(params::Dict)

Returns the name of the mesh file from the parameters

# Arguments
- `params::Dict`: The parameters
# Returns
- `String`: The name of the finite element topology file
"""
function get_external_topology_name(params::Dict)
    check = haskey(params["Discretization"], "Input External Topology")
    if !check
        return nothing
    end
    check = haskey(params["Discretization"]["Input External Topology"], "File")
    if !check
        @error "Input External Topology is defined without a file where to find it."
        return nothing
    end
    return params["Discretization"]["Input External Topology"]["File"]
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
    get_header(filename::Union{String,AbstractString})

Returns the header line and the header.

# Arguments
- `filename::Union{String,AbstractString}`: The filename of the file.
# Returns
- `header_line::Int`: The header line.
- `header::Vector{String}`: The header.
"""
function get_header(filename::Union{String,AbstractString})
    file = open(filename, "r")
    header_line = 0
    for line in eachline(file)#
        header_line += 1
        if contains(line, "header:")
            close(file)
            return header_line, convert(Vector{String}, split(line)[2:end])
        end
    end
    @error "No header exists in $filename. Please insert 'header: global_id' above the first node"
    return
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
    type = get(params["Discretization"], "Type", "Text File")
    if type == "Exodus"
        exo = ExodusDatabase(joinpath(path, get_mesh_name(params)), "r")
        nset_names = read_names(exo, NodeSet)
        conn = collect_element_connectivities(exo)
        nset_nodes = []
        for (id, entry) in enumerate(nset_names)
            nset_nodes = Vector{Int64}(read_set(exo, NodeSet, id).nodes)
            if length(entry) == 0
                nsets["Set-"*string(id)] = findall(row -> all(val -> any(val .== nset_nodes), row), conn)
            else
                nsets[entry] = findall(row -> all(val -> any(val .== nset_nodes), row), conn)
            end
            # end
        end
        conn = nothing
        nset_nodes = nothing
        @info "Found $(length(nsets)) node sets"
        close(exo)
        return nsets
    end
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