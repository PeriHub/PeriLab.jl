# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using LinearAlgebra
using AbaqusReader
using DataFrames
using CSV

"""
    extrude_surface_mesh(mesh::DataFrame)

extrude the mesh at the surface of the block

# Arguments
- `mesh::DataFrame`: The input mesh data represented as a DataFrame.
- `params::Dict`: The input parameters.
"""
function extrude_surface_mesh(mesh::DataFrame, direction, step, number)

    # Finding min and max values for each dimension
    min_x, max_x = extrema(mesh.x)
    min_y, max_y = extrema(mesh.y)
    min_z = 0.0
    max_z = 0.0
    if "z" in names(mesh)
        min_z, max_z = extrema(mesh.z)
    end

    if direction == "X"
        coord_min = min_x
        coord_max = max_x
        row_min = min_y
        row_max = max_y
    elseif direction == "Y"
        coord_min = min_y
        coord_max = max_y
        row_min = min_x
        row_max = max_x
    end

    block_id = maximum(mesh.block_id) + 1
    volume = step * step

    id = 0

    node_sets = Dict("Extruded_1" => [], "Extruded_2" => [])

    for i in coord_max+step:step:coord_max+step*number, j in row_min:step:row_max+step, k in min_z:step:max_z+step
        if direction == "X"
            push!(mesh, (x=i, y=j, z=k, volume=volume, block_id=block_id))
        elseif direction == "Y"
            push!(mesh, (x=j, y=i, z=k, volume=volume, block_id=block_id))
        end
        append!(node_sets["Extruded_1"], [Int64(id)])
        id += 1
    end

    block_id += 1

    for i in coord_min-step:-step:coord_min-step*number, j in row_min:step:row_max+step, k in min_z:step:max_z+step
        if direction == "X"
            push!(mesh, (x=i, y=j, z=k, volume=volume, block_id=block_id))
        elseif direction == "Y"
            push!(mesh, (x=j, y=i, z=k, volume=volume, block_id=block_id))
        end
        append!(node_sets["Extruded_2"], [Int64(id)])
        id += 1
    end
    return mesh, node_sets
end

"""
    calculate_volume(element_type::String, vertices::Vector{Vector{Float64}})

Calculate the volume of a element.

# Arguments
- `element_type`: The element type of the element.
- `vertices`: The vertices of the element.
# Returns
- `volume`: The volume of the element.
"""
function calculate_volume(element_type::String, vertices::Vector{Vector{Float64}})
    if element_type == "Quad4"
        return area_of_polygon(vertices)
    elseif element_type == "Tet4"
        return tetrahedron_volume(vertices)
    elseif element_type == "Wedge6"
        return wedge6_volume(vertices)
    elseif element_type == "Hex8"
        return hex8_volume(vertices)
    else
        @error "Element type $element_type currently not supported"
        return nothing
    end
end

"""
    tetrahedron_volume(tet_vertices)

Calculate the volume of a tetrahedron.

# Arguments
- `tet_vertices`: The vertices of the tetrahedron.
# Returns
- `volume`: The volume of the tetrahedron.
"""
function tetrahedron_volume(tet_vertices::Vector{Vector{Float64}})
    mat = hcat(hcat(tet_vertices...)', ones(4))  # Augmenting matrix with ones in the fourth column
    volume = abs(det(mat) / 6)   # Using det function to calculate determinant
    return volume
end

"""
hex8_volume(hex_vertices)

Calculate the volume of a hex.

# Arguments
- `hex_vertices`: The vertices of the wedge.
# Returns
- `volume`: The volume of the wedge.
"""
function hex8_volume(hex_vertices::Vector{Vector{Float64}})
    tets = [
        [hex_vertices[1], hex_vertices[2], hex_vertices[4], hex_vertices[5]],
        [hex_vertices[2], hex_vertices[3], hex_vertices[4], hex_vertices[7]],
        [hex_vertices[2], hex_vertices[5], hex_vertices[6], hex_vertices[7]],
        [hex_vertices[4], hex_vertices[5], hex_vertices[7], hex_vertices[8]],
        [hex_vertices[2], hex_vertices[4], hex_vertices[5], hex_vertices[7]]
    ]

    volumes = []
    for tet in tets
        push!(volumes, tetrahedron_volume(tet))
    end

    return sum(volumes)
end

"""
    wedge6_volume(wedge_vertices)

Calculate the volume of a wedge.

# Arguments
- `wedge_vertices`: The vertices of the wedge.
# Returns
- `volume`: The volume of the wedge.
"""
function wedge6_volume(wedge_vertices::Vector{Vector{Float64}})
    tets = [
        [wedge_vertices[1], wedge_vertices[2], wedge_vertices[3], wedge_vertices[4]],
        [wedge_vertices[2], wedge_vertices[3], wedge_vertices[4], wedge_vertices[5]],
        [wedge_vertices[3], wedge_vertices[4], wedge_vertices[5], wedge_vertices[6]]
    ]

    volumes = []
    for tet in tets
        push!(volumes, tetrahedron_volume(tet))
    end

    return sum(volumes)
end

"""
    area_of_polygon(vertices)

Calculate the area of a polygon.

# Arguments
- `vertices`: The vertices of the polygon.
# Returns
- `area`: The area of the polygon.
"""
function area_of_polygon(vertices)
    n = length(vertices)
    area = 0.0

    for i in 1:n
        j = mod(i, n) + 1
        area += (vertices[i][1] + vertices[j][1]) * (vertices[i][2] - vertices[j][2])
    end

    return abs(area) / 2.0
end

function read(filename)
    mesh = abaqus_read_mesh(filename; verbose=false)

    nodes = mesh["nodes"]
    elements = mesh["elements"]
    element_sets = mesh["element_sets"]
    element_types = mesh["element_types"]

    dof = 2

    nodes_vector = collect(values(nodes))
    node_z = [node[end] for node in nodes_vector]
    if length(unique(node_z)) > 2
        dof = 3
    end
    @info "Abaqus mesh with $dof DOF"

    mesh_df = ifelse(dof == 2,
        DataFrame(x=[], y=[], volume=[], block_id=Int[]),
        DataFrame(x=[], y=[], z=[], volume=[], block_id=Int[])
    )

    id = 1
    block_id = 1
    element_written = []
    nsets = Dict{String,Any}()

    # sort element_sets by length
    # element_sets_keys = sort(collect(keys(element_sets)), by=x -> length(element_sets[x]), rev=true)
    element_sets_keys = collect(keys(element_sets))
    for key in element_sets_keys
        ns_nodes = Int64[]
        for (i, element_id) in enumerate(element_sets[key])
            if element_id in element_written
                push!(ns_nodes, findfirst(x -> x == element_id, element_written))
                continue
            end
            push!(ns_nodes, id)
            node_ids = elements[element_id]
            element_type = element_types[element_id]
            vertices = [nodes[node_id] for node_id in node_ids]
            volume = calculate_volume(string(element_type), vertices)
            center = sum(vertices) / size(vertices)[1]
            if dof == 2
                push!(mesh_df, (x=center[1], y=center[2], volume=volume, block_id=block_id))
            else
                push!(mesh_df, (x=center[1], y=center[2], z=center[3], volume=volume, block_id=block_id))
            end
            push!(element_written, element_id)
            id += 1
        end
        nsets[key] = ns_nodes
        block_id += 1
    end
    @info "Found $(maximum(mesh_df.block_id)) block(s)"
    @info "Found $(length(nsets)) node sets"
    @info "NodeSets: $element_sets_keys"

    nnodes = size(mesh_df, 1) + 1

    mesh_df, surface_ns = extrude_surface_mesh(mesh_df, "X", 3, 2)

    for (key, values) in surface_ns
        nsets[key] = Vector{Int64}(values .+ nnodes)
    end

    txt_file = replace(filename, ".inp" => ".txt")
    write(txt_file, "header: x y volume block_id\n")
    CSV.write(txt_file, mesh_df; delim=' ', append=true)
    for (key, values) in nsets
        txt_file = replace("ns_" * filename, ".inp" => "_" * key * ".txt")
        write(txt_file, "header: global_id\n")
        for value in values
            write(txt_file, "$value\n")
        end
    end

    mesh = nothing
    nodes = nothing
    elements = nothing
    element_sets = nothing

end

read("INPUTFILE.inp")