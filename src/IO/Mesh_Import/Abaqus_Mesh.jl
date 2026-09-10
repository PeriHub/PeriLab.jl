# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
"""
    Abaqus_Mesh

Abaqus mesh importer. Reads an Abaqus INP file and converts the element
connectivities into the peridynamic mesh representation (center of each
element as a node, element volume, block id).
"""
module Abaqus_Mesh

using ...PeriLabExceptions: @abort
using DataFrames
using AbaqusReader
using .Mesh_Volume: calculate_volume

export mesh_import_name
export read_mesh

"""
    mesh_import_name()

Gives the mesh importer name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: The name of the mesh importer.
"""
function mesh_import_name()
    return "Abaqus"
end

"""
    read_mesh(params::Dict, filename::String)

Reads an Abaqus mesh file and returns the mesh data as a DataFrame along with
the node sets found in the boundary conditions.

# Arguments
- `params::Dict`: The parameters.
- `filename::String`: The path to the Abaqus INP file.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
- `nsets::Dict{String,Vector{Int64}}`: The node sets.
"""
function read_mesh(params::Dict, filename::String)
    mesh = abaqus_read_mesh(filename; verbose = false)

    nodes = mesh["nodes"]
    elements = mesh["elements"]
    element_sets = mesh["element_sets"]
    @assert length(element_sets) > 0
    element_types = mesh["element_types"]

    dof = 2
    nodes_vector = collect(values(nodes))
    if size(nodes_vector[1])[1] == 3
        dof = 3
    end
    @info "Abaqus mesh with $dof DOF"

    num_elements = length(elements)
    mesh_df = ifelse(dof == 2,
                     DataFrame(x = Array{Float64,1}(undef, num_elements),
                               y = Array{Float64,1}(undef, num_elements),
                               volume = Array{Float64,1}(undef, num_elements),
                               block_id = Array{Int64,1}(undef, num_elements)),
                     DataFrame(x = Array{Float64,1}(undef, num_elements),
                               y = Array{Float64,1}(undef, num_elements),
                               z = Array{Float64,1}(undef, num_elements),
                               volume = Array{Float64,1}(undef, num_elements),
                               block_id = Array{Int64,1}(undef, num_elements)))

    id = 1
    block_id = 1
    element_written = []
    nsets = Dict{String,Vector{Int64}}()

    nset_names = []

    for boundary_condtion in keys(params["Boundary Conditions"])
        if haskey(params["Boundary Conditions"][boundary_condtion], "Node Set")
            push!(nset_names,
                  params["Boundary Conditions"][boundary_condtion]["Node Set"])
        end
    end
    nset_names = unique(nset_names)

    element_sets_keys = collect(keys(element_sets))
    for nset in nset_names
        if nset in element_sets_keys
            deleteat!(element_sets_keys, findfirst(x -> x == nset, element_sets_keys))
            push!(element_sets_keys, nset)
        end
    end
    block_names = []
    for key in element_sets_keys
        element_set = element_sets[key]
        ns_nodes = Array{Int64,1}(undef, length(element_set))
        nset_only = true
        for (jID, element_id) in enumerate(element_set)
            if element_id in element_written
                if key in nset_names
                    ns_nodes[jID] = findfirst(x -> x == element_id, element_written)
                end
                continue
            end
            nset_only = false
            ns_nodes[jID] = id
            node_ids = elements[element_id]
            element_type = element_types[element_id]
            vertices = [nodes[node_id] for node_id in node_ids]
            volume = calculate_volume(string(element_type), vertices)
            center = sum(vertices) / length(vertices)
            if dof == 2
                mesh_df[id, :] = [center[1], center[2], volume, block_id]
            else
                mesh_df[id, :] = [center[1], center[2], center[3], volume, block_id]
            end
            push!(element_written, element_id)
            id += 1
        end
        if key in nset_names
            nsets[key] = ns_nodes
        end
        if !nset_only
            block_id += 1
            push!(block_names, key)
        end
    end
    @info "Found $(block_id-1) block(s)"
    @info "Blocks: $block_names"
    @info "Found $(length(nsets)) node set(s)"
    if length(nsets) > 0
        @info "NodeSets: $(keys(nsets))"
    end

    mesh = nothing
    nodes = nothing
    elements = nothing
    element_sets = nothing

    return mesh_df, nsets
end

end # module Abaqus_Mesh
