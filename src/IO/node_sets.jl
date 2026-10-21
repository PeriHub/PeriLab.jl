# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using CSV
using Exodus
using ..InputDeck: DiscretizationParams

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
    @warn "No header exists in $filename. Please insert 'header: global_id' above the first node"
end

"""
    external_topology_file(d, path)

File name of the external element topology, or `nothing`. Aborts if the file
does not exist in `path`.
"""
function external_topology_file(d::DiscretizationParams, path::String)
    topology = d.input_external_topology
    topology === nothing && return nothing
    filename = joinpath(path, topology.file)
    if !isfile(filename)
        @abort "External topology file: ''$filename'' does not exist"
        return
    end
    return topology.file
end

"""
    read_node_sets(d, path, mesh_df)

Node sets of the discretization: from an Exodus mesh, or from `Node Sets`
entries (node id, id list, range `a:b`, coordinate expression in `x`/`y`/`z`,
`All`, or a node set file).
"""
function read_node_sets(d::DiscretizationParams, path::String, mesh_df::DataFrame)
    nsets = Dict{String,Vector{Int64}}()
    if d.type == "Exodus"
        exo = ExodusDatabase(joinpath(path, d.input_mesh_file), "r")
        nset_names = read_names(exo, NodeSet)
        conn = collect_element_connectivities(exo)
        for (id, entry) in enumerate(nset_names)
            nset_nodes = Vector{Int64}(read_set(exo, NodeSet, id).nodes)
            name = length(entry) == 0 ? "Set-" * string(id) : entry
            nsets[name] = findall(row -> all(val -> any(val .== nset_nodes), row), conn)
        end
        @info "Found $(length(nsets)) node sets"
        close(exo)
        return nsets
    end
    for (entry, value) in d.node_sets
        if value isa Int64
            nsets[entry] = [value]
        elseif occursin(".txt", value)
            if isnothing(get_header(joinpath(path, value)))
                @warn "Node set file " * value *
                      " is not correctly specified. Please check the examples. The node set is excluded."
                continue
            end
            header_line, header = get_header(joinpath(path, value))
            nodes = CSV.read(joinpath(path, value), DataFrame; delim = " ", header = false,
                             skipto = header_line + 1,)
            if size(nodes) == (0, 0)
                @abort "Node set file is empty " * value * ". The node set is excluded."
                return
            end
            nsets[entry] = nodes.Column1
        elseif occursin(":", value)
            nsets[entry] = collect(eval(Meta.parse(value)))
        elseif occursin("x", value) || occursin("y", value) || occursin("z", value)
            nodes = []
            for id in 1:size(mesh_df, 1)
                global x = mesh_df[!, "x"][id]
                global y = mesh_df[!, "y"][id]
                if occursin("z", value)
                    global z = mesh_df[!, "z"][id]
                end
                try
                    if eval(Meta.parse(value))
                        push!(nodes, id)
                    end
                catch UndefVarError
                    @abort "Failed to eval nodeset value: '$(value)', $UndefVarError"
                    return
                end
            end
            nsets[entry] = nodes
        elseif value == "All"
            nsets[entry] = collect(1:size(mesh_df, 1))
        else
            nsets[entry] = parse.(Int, split(value))
        end
    end
    return nsets
end
