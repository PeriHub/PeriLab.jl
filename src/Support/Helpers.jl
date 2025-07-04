# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Helpers

using PointNeighbors: GridNeighborhoodSearch, initialize_grid!, foreach_neighbor
using Meshes
using Tensors
using Dierckx
using ProgressBars
using LinearAlgebra
using StaticArrays
using LoopVectorization
using Unitful
using CDDLib, Polyhedra

export qdim
export check_inf_or_nan
export compute_geometry
export find_active_nodes
export get_active_update_nodes
export find_indices
export find_inverse_bond_id
export find_files_with_ending
export get_block_nodes
export matrix_style
export get_fourth_order
export interpolation
export interpol_data
export progress_bar
export invert
export compute_distance_and_normals
export rotate
export sub_in_place!
export add_in_place!
export div_in_place!
export fill_in_place!
export fastdot
export fast_mul!
export mat_mul!
export get_mapping
export mat_mul_transpose_mat!
export get_ring
export get_hexagon
export nearest_point_id
export get_shared_horizon

function get_shared_horizon(datamanager, id)
    horizon = datamanager.get_field("Shared Horizon")
    return horizon[id]
end

"""
     matrix_to_vector(mat::Union{Matrix{Float64},Matrix{Int64}})

Transforming a matrix representation in a Vector{Vector} representation.

# Arguments
- `mat::Union{Matrix{Float64},Matrix{Int64}}`: Points which form the polyhedron.
# Returns
- ``: transformed data
"""
function matrix_to_vector(mat::Union{Matrix{Float64},Matrix{Int64}})
    return [vec(mat[i, :]) for i in eachindex(axes(mat, 1))]
end

"""
    compute_geometry(points::Union{Matrix{Float64},Matrix{Int64}})

Returns a polyhedron object in 2D or 3D. To do so the matrix form of storing the geometry is transfered in a Vector{Vector}. This form is than transformed in a V-representation of the Polyhedra.jl package.

# Arguments
- `points::Union{Matrix{Float64},Matrix{Int64}}`: Points which form the polyhedron.
# Returns
- ``: polyhedron
"""

function compute_geometry(points::Union{Matrix{Float64},Matrix{Int64}})
    return polyhedron(vrep(matrix_to_vector(points)))
end

"""
    compute_surface_nodes_and_connections(points::Union{Matrix{Float64},Matrix{Int64}},
                                               poly, free_surfaces::Vector{Int64})

Computes the points which are connected free surfaces (()).
This function is used for contact search purposes. The free surface nodes are used to compute the nearest neighbors. The connections and underlying points are needed for the contact algorithm. They are a subset of the surface points to create the polyhedron and the only ones which can be in contact.

# Arguments
- `points::Union{Matrix{Float64},Matrix{Int64}}`: Points which form the polyhedron.
- `poly`: Polyhedron object.
- `free_surfaces::Vector{Int64}`: List of the free surfaces of the polyhedron.

# Returns
- `connections`: Tthe connections to the free surfaces. There can be more surface points than connections.
"""
function compute_free_surface_nodes(points::Union{Matrix{Float64},Matrix{Int64}},
                                    ids::Vector{Int64},
                                    free_surfaces::Vector{Int64})
    poly = compute_geometry(points[ids, :])
    normals, offset = get_surface_information(poly)
    free_nodes = Vector{Int64}([])
    for pID in ids
        for id in free_surfaces
            if length(offset) < id #TODO: check this condition
                @error "Surface ID $id is not defined"
            end
            if isapprox(dot(normals[id, :], points[pID, :]) - offset[id], 0; atol = 1e-6)
                # connections only to the free surfaces.
                append!(free_nodes, pID)
            end
        end
    end
    # not possible; deformation leads to a switch in surface numbers
    return sort(free_nodes)
end

# taken from peridigm
function compute_distance_and_normals(p1, p2)
    distance = norm(p2 - p1)
    distance == 0 ? distance + eps() :
    distance
    return distance, (p2 - p1) ./ distance
end

function get_surface_normals(poly)
    return MixedMatHRep(hrep(poly)).A
end

function get_surface_offset(poly)
    return MixedMatHRep(hrep(poly)).b
end

function get_surface_information(poly)
    normals = get_surface_normals(poly)
    b = get_surface_offset(poly)
    for i in eachindex(b)
        b[i] /= norm(normals[i, :])
        normals[i, :] ./= norm(normals[i, :])
    end
    return normals, b
end

function check_neighbor_position(poly, points, nlist::Vector{Int64}, msg::Bool)
    for nID in nlist
        if !point_is_inside(points[nID, :], poly)
            if msg
                msg = false
                @warn "Make sure that your contact block is large enough. If it is surrounded by non contact blocks some of the points near the edgdes are ignored."
            end
            return false
        end
    end
    return true
end

function point_is_inside(point, poly)
    return point in poly
end

function nearest_point_id(p, points)
    return argmin(norm.(eachrow(points) .- Ref(p)))
end

"""
    remove_ids(dict::Dict{Int64,Int64}, numbers::Vector{Int64})

Remove multiple keys in a Vector from a Dict.

# Arguments
- `dict::Dict{Int64,Int64}`: Dictionary with ids
- `numbers::Vector{Int64}`: Ids to remove
# Returns
updated dict
"""

function remove_ids(dict::Dict{Int64,Int64}, numbers::Vector{Int64})
    for num in numbers
        delete!(dict, num)
    end
end

function remove_ids(vec::Vector{Int64}, numbers::Vector{Int64})
    filter!(x -> !(x in numbers), vec)
end

"""
    get_block_nodes(block_ids, nnodes)

Returns a dictionary mapping block IDs to collections of nodes.

# Arguments
- `block_ids::Vector{Int64}`: A vector of block IDs
- `nnodes::Int64`: The number of nodes
# Returns
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes
"""
function get_block_nodes(block_ids, nnodes)
    block_nodes = Dict{Int64,Vector{Int64}}()
    for i in unique(block_ids[1:nnodes])
        block_nodes[i] = find_indices(block_ids[1:nnodes], i)
    end
    return block_nodes
end

function get_nearest_neighbors(nodes,
                               dof::Int64,
                               system_coordinates,
                               neighbor_coordinates,
                               radius::Union{Int64,Float64,Vector{Float64},Vector{Int64}},
                               neighborList,
                               diffent_lists = false)
    nhs = GridNeighborhoodSearch{dof}(search_radius = maximum(radius),
                                      n_points = length(nodes))
    initialize_grid!(nhs, system_coordinates')
    initialize_grid!(nhs, neighbor_coordinates')
    list_empty = true

    for iID in nodes
        neighbors = []
        foreach_neighbor(system_coordinates',# System coordinates -> must be transpose for the datamanager definition
                         neighbor_coordinates',# Potential neighbor coordinates -> must be transpose for the datamanager definition
                         nhs,
                         iID,
                         search_radius = radius isa AbstractVector ? radius[iID] : radius) do i,
                                                                                              j,
                                                                                              _,
                                                                                              L
            if i != j || diffent_lists
                push!(neighbors, j)
                list_empty = false
            end
        end
        neighborList[iID] = neighbors
    end
    return neighborList, list_empty
end

function find_point_in_polygon(point, coor, topo)
    return Point(point[1], point[2]) in get_ring(coor, topo)
end

function find_point_in_hexagon(point, coor, topo)
    return Point(point[1], point[2], point[3]) in get_hexagon(coor, topo)
end

"""
    get_hexagon(coor::Union{Matrix{Float64},Matrix{Int64}}, topo::Vector{Int64})

# Description
This function gives the hexahedron model back to compute centroids of this surface and check if a point lies inside. The el_topology has to be transformed from PeriLab notation, to avoid overlapping.

# Arguments
- `coor::Union{Matrix{Float64},Matrix{Int64}}`: Coordinates.
- `topo::Vector{Int64}`: el_topology of an element.
# Returns
- `Hexahedron`: eight point hexahedron.

"""
function get_hexagon(coor::Union{Matrix{Float64},Matrix{Int64}}, topo::Vector{Int64})
    return Hexahedron(Point(coor[topo[1], 1], coor[topo[1], 2], coor[topo[1], 3]),
                      Point(coor[topo[2], 1], coor[topo[2], 2], coor[topo[2], 3]),
                      Point(coor[topo[4], 1], coor[topo[4], 2], coor[topo[4], 3]),
                      Point(coor[topo[3], 1], coor[topo[3], 2], coor[topo[3], 3]),
                      Point(coor[topo[5], 1], coor[topo[5], 2], coor[topo[5], 3]),
                      Point(coor[topo[6], 1], coor[topo[6], 2], coor[topo[6], 3]),
                      Point(coor[topo[8], 1], coor[topo[8], 2], coor[topo[8], 3]),
                      Point(coor[topo[7], 1], coor[topo[7], 2], coor[topo[7], 3]))
end

"""
    get_ring(coor::Union{Matrix{Float64},Matrix{Int64}}, topo::Vector{Int64})

# Description
This function gives the ring model back to compute centroids of this surface and check if a point lies inside. The el_topology has to be transformed from PeriLab notation, to avoid overlapping.

# Arguments
- `coor::Union{Matrix{Float64},Matrix{Int64}}`: Coordinates.
- `topo::Vector{Int64}`: el_topology of an element.
# Returns
- `Ring`: four point closed surface.

"""
function get_ring(coor::Union{Matrix{Float64},Matrix{Int64}}, topo::Vector{Int64})
    return Ring(Point(coor[topo[1], 1], coor[topo[1], 2]),
                Point(coor[topo[2], 1], coor[topo[2], 2]),
                Point(coor[topo[4], 1], coor[topo[4], 2]),
                Point(coor[topo[3], 1], coor[topo[3], 2]))
end

function get_centroid(el, dof)
    c = centroid(el)
    if dof == 2
        return ustrip(c.coords.x), ustrip(c.coords.y)
    else
        return ustrip(c.coords.x), ustrip(c.coords.y), ustrip(c.coords.z)
    end
end
"""
    create_centroid_and_search_radius(coor::Union{Matrix{Float64},Matrix{Int64}},
                                      el_topology::Matrix{Int64},
                                      dof::Int64,
                                      fu)

Computes the centroid and Contact Radius for each element in a given mesh.

# Arguments
- `coor::Union{Matrix{Float64}, Matrix{Int64}}`:   A matrix of size `(N x dof)`, where each row represents the coordinates of a point in the space.
- `el_topology::Matrix{Int64}`:   A matrix of size `(M x num_nodes)`, where each row contains the indices of points that define an element.
- `dof::Int64`:   The number of spatial dimensions (degrees of freedom), e.g., `2` for 2D and `3` for 3D.
- `fu`:   A function of four point 2D surface or eight point 3D volume.

# Returns
- `el_centroid::Matrix{Float64}`:
  A matrix of size `(M × dof)`, where each row represents the centroid of an element.
- `search_radius::Vector{Float64}`:
  A vector of size `M`, where each entry is the maximum distance from the element centroid to any of its nodes.
"""
function create_centroid_and_search_radius(coor::Union{Matrix{Float64},Matrix{Int64}},
                                           el_topology::Matrix{Int64},
                                           dof::Int64,
                                           fu)
    el_centroid = zeros(Float64, size(el_topology)[1], dof)
    search_radius = zeros(Float64, size(el_topology)[1])
    for iel in eachindex(el_topology[:, 1])
        topo = el_topology[iel, :]
        el_centroid[iel, :] .= get_centroid(fu(coor, topo), dof)
        search_radius[iel] = maximum([norm(el_centroid[iel, :] .- p) for p in coor[topo, :]])
    end
    return el_centroid, search_radius
end

function find_point_in_element(el_topology, near_points, coor, fu, coupling_dict)
    @inbounds @fastmath for iel in axes(el_topology, 1)
        topo = el_topology[iel, :]
        for point in near_points[iel]
            if haskey(coupling_dict, point)
                # only of elment per point for coupling
                continue
            end
            if fu(coor[point, :], coor, topo)
                coupling_dict[point] = iel
            end
        end
    end
    return coupling_dict
end

function get_mapping(dof::Int64)
    if dof == 2
        return (1, 1), (2, 2), (2, 1)
    elseif dof == 3
        return (1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2)
    else
        @error "$dof is no valid mapping option."
        return nothing
    end
end

function matrix_diff!(s3::AbstractArray{Float64,3}, nodes::AbstractArray{Int64},
                      s2::AbstractArray{Float64,3},
                      s1::AbstractArray{Float64,3})
    @views @inbounds for iID in nodes
        @inbounds @fastmath @views for m in axes(s1[iID, :, :], 1),
                                       n in axes(s1[iID, :, :], 2)
            s3[iID, m, n] = s2[iID, m, n] - s1[iID, m, n]
        end
    end
end

function fast_mul!(stress_NP1::AbstractMatrix{Float64},
                   C::AbstractArray{Float64},
                   strain_increment::AbstractMatrix{Float64},
                   stress_N::AbstractMatrix{Float64},
                   mapping::NTuple{N,Tuple{Int,Int}}) where {N}
    @inbounds @fastmath for m in 1:N
        i, j = mapping[m]
        sNP1 = stress_N[i, j]

        for k in 1:N
            ki, kj = mapping[k]
            sNP1 += C[m, k] * strain_increment[ki, kj]
        end

        stress_NP1[i, j] = sNP1
        # Only set symmetric component if i != j
        if i != j
            stress_NP1[j, i] = sNP1
        end
    end

    return stress_NP1
end

function mat_mul!(C::AbstractMatrix{T}, A::AbstractMatrix{T},
                  B::AbstractMatrix{T}) where {T<:Float64}
    @inbounds @fastmath for m in axes(A, 1), n in axes(B, 2)
        Cmn = zero(T)
        for k in axes(A, 2)
            Cmn += A[m, k] * B[k, n]
        end
        C[m, n] = Cmn
    end
end
function mat_mul_transpose_mat!(C::AbstractMatrix{T}, A::AbstractMatrix{T},
                                B::AbstractMatrix{T}) where {T<:Float64}
    @inbounds @fastmath for m in axes(A, 1), n in axes(B, 2)
        Cmn = zero(T)
        @inbounds @fastmath for k in axes(A, 2)
            Cmn += A[k, m] * B[k, n]
        end
        C[m, n] = 0.5 * Cmn
    end
end
function add_in_place!(C::AbstractMatrix{T},
                       A::Vector{Vector{T}},
                       B::Vector{Vector{T}},
                       factor = 1) where {T<:Union{Int64,Float64,Bool}}
    m = length(A)
    n = length(A[1])

    for i in 1:n
        for j in 1:n
            for k in 1:m
                C[i, j] += A[k][i] * B[k][j] * factor
            end
        end
    end

    return C
end
function get_MMatrix(len::Int64)::MMatrix
    if len == 4
        return MMatrix{2,2}(zeros(Float64, 2, 2))
    elseif len == 9
        return MMatrix{3,3}(zeros(Float64, 3, 3))
    elseif len == 36
        return MMatrix{6,6}(zeros(Float64, 6, 6))
    else
        @error "MMatrix length $len not pre-allocated. Please add your size to helper.jl in get_MMatrix."
        return MMatrix{1,1}(zeros(Float64, 1, 1))
    end
end
function mul_in_place!(C::Vector{T}, A::Vector{T},
                       B::Vector{Bool}) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        C[i] = A[i] * B[i]
    end
end
function mul_in_place!(C::Vector{Vector{T}},
                       A::Vector{Vector{T}},
                       B::Vector{T}) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        @inbounds for j in eachindex(A[i])
            C[i][j] = A[i][j] * B[i]
        end
    end
end
function mul_in_place!(D::Vector{Vector{T}},
                       A::Vector{Vector{T}},
                       B::Vector{Vector{T}},
                       C::T = 1.0) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(D) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        @inbounds for j in eachindex(A[i])
            D[i][j] = A[i][j] * B[i][j] * C
        end
    end
end
function sub_in_place!(C::Vector{Vector{T}},
                       A::Vector{Vector{T}},
                       B::Vector{Vector{T}}) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        C[i] .= A[i] .- B[i]
    end
end
function sub_in_place!(C::Vector{Vector{Vector{T}}},
                       A::Vector{Vector{Vector{T}}},
                       B::Vector{Vector{Vector{T}}}) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        @inbounds for j in eachindex(A[i])
            C[i][j] .= A[i][j] .- B[i][j]
        end
    end
end
function add_in_place!(C::Vector{Vector{T}},
                       A::Vector{Vector{T}},
                       B::Vector{Vector{T}}) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        C[i] .= A[i] .+ B[i]
    end
end
function div_in_place!(C::Vector{T}, A::Vector{T}, B::T,
                       absolute = false) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(C) == length(A)

    @inbounds for i in eachindex(A)
        if absolute
            C[i] = abs(A[i]) / B
        else
            C[i] = A[i] / B
        end
    end
end
function div_in_place!(C::Vector{Vector{T}},
                       A::Vector{Vector{T}},
                       B::Vector{Vector{T}}) where {T<:Union{Int64,Float64,Bool}}
    # Check if dimensions match
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        C[i] .= A[i] ./ B[i]
    end
end
function fastdot(a::AbstractArray, b::AbstractArray, absolute = false)
    c = zero(eltype(a))
    @assert length(a) == length(b)
    @inbounds @simd for i in eachindex(a, b)
        if absolute
            c += abs(a[i]) * abs(b[i])
        else
            c += a[i] * b[i]
        end
    end
    c
end
function fill_in_place!(A::Vector{Vector{T}},
                        value::T,
                        active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool}}
    @inbounds for i in eachindex(A)
        if active[i]
            A[i] .= value
        end
    end
end

function fill_in_place!(A::Vector{Array{T,N}},
                        value::T,
                        active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool},N}
    @inbounds for i in eachindex(A)
        if active[i]
            A[i] .= value
        end
    end
end
function fill_in_place!(A::Vector{Vector{Vector{T}}},
                        value::T,
                        active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool}}
    @inbounds for i in eachindex(A)
        if active[i]
            @inbounds for j in eachindex(A[i])
                A[i][j] .= value
            end
        end
    end
end
function fill_in_place!(A::Vector{Vector{Array{T,N}}},
                        value::T,
                        active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool},N}
    @inbounds for i in eachindex(A)
        if active[i]
            @inbounds for j in eachindex(A[i])
                A[i][j] .= value
            end
        end
    end
end
# function fourth_order_times_second_order_tensor(C, s1, s2, s3, dof)

#     @inbounds @simd for i = 1:dof
#         for j = 1:dof
#             s1[i, j] = s3[i, j] + fastdot(C[i, j, :, :], s2[:, :])
#         end
#     end
#     return s1
# end

"""
    qdim(order::Int64, dof::Int64)

Calculate the number of terms in a polynomial expansion up to a specified accuracy order. Simplied first complex loop in Peridigm correspondence::computeLagrangianGradientWeights.
In the unit test this values where tested.

# Arguments
- `order::Int64`: The accuracy order of the polynomial expansion.

# Returns
- `Int64`: The total number of terms in the polynomial expansion.

# Description
This function calculates the number of terms in a polynomial expansion up to the specified accuracy order
using an analytical formula derived from combinatorial considerations. The function iterates over each order
from 1 to the specified `order` and calculates the sum of binomial coefficients according to the formula:
qdim(order) = Σ(i=1 to order) [(i+2)! / (2! * i!)]

"""
function qdim(order::Int64, dof::Int64)
    if order < 1
        @error "Accuracy order must be greater than zero."
        return nothing
    end
    return sum(binomial(i + dof - 1, dof - 1) for i in 1:order)
end

"""
    find_indices(vector, what)

Returns the indices of `vector` that are equal to `what`.

# Arguments
- `vector::Vector`: The vector to search in.
- `what`: The value to search for.
# Returns
- `indices::Vector`: The indices of `vector` that are equal to `what`.
"""
function find_indices(vector::Vector{T}, what::T) where {T<:Union{Int64,Float64,Bool}}
    return findall(item -> item == what, vector)
end

"""
    get_active_update_nodes(active::Vector{Bool}, update_list::Vector{Bool}, nodes::Vector{Int64}, index::Vector{Int64})

Returns the active nodes and the update nodes.

# Arguments
- `active::Vector{Bool}`: The active vector.
- `update_list::Vector{Bool}`: The update vector.
- `nodes::Vector{Int64}`: The vector of nodes.
- `index::Vector{Int64}`: Pre allocated Vector.
# Returns
- `update_nodes::Vector{Int64}`: The nodes of `update` that are true.
"""
function get_active_update_nodes(active, update, nodes, index)
    count::Int64 = 0
    for node in nodes
        if active[node] && update[node]
            count += 1
            index[count] = node
        end
    end
    return view(index, 1:count)
end

function find_active_nodes(active,
                           active_nodes::Union{Vector{Int64},SubArray},
                           nodes,
                           false_or_true::Bool = true)
    count::Int64 = 0
    for node in nodes
        if active[node] == false_or_true
            count += 1
            active_nodes[count] = node
        end
    end
    return view(active_nodes, 1:count)
end

"""
    find_files_with_ending(folder_path::AbstractString, file_ending::AbstractString)

Returns a list of files in `folder_path` that end with `file_ending`.

# Arguments
- `folder_path::AbstractString`: The path to the folder.
- `file_ending::AbstractString`: The ending of the files.
# Returns
- `file_list::Vector{String}`: The list of files that end with `file_ending`.
"""
function find_files_with_ending(folder_path::AbstractString, file_ending::AbstractString)
    file_list = filter(x -> isfile(joinpath(folder_path, x)) && endswith(x, file_ending),
                       readdir(folder_path))
    return file_list
end

"""
    check_inf_or_nan(array, msg)

Checks if the sum of the array is finite. If not, an error is raised.

# Arguments
- `array`: The array to check.
- `msg`: The error message to raise.
# Returns
- `Bool`: `true` if the sum of the array is finite, `false` otherwise.
"""
function check_inf_or_nan(array::AbstractArray{T},
                          msg::String) where {T<:Union{Int64,Float64}}
    if isnan(sum(array))
        @error "Field ''$msg'' has NaN elements."
        return true
    end
    if !isfinite(sum(array))
        @error "Field ''$msg'' is infinite."
        return true
    end
    return false
end
function check_inf_or_nan(array::T,
                          msg::String) where {T<:Union{Int64,Float64}}
    if isnan(sum(array))
        @error "Field ''$msg'' has NaN elements."
        return true
    end
    if !isfinite(sum(array))
        @error "Field ''$msg'' is infinite."
        return true
    end
    return false
end
"""
    matrix_style(A)

Include a scalar or an array and reshape it to style needed for LinearAlgebra package

# Arguments
- `A`: The array or scalar to reshape
# Returns
- `Array`: The reshaped array
"""
function matrix_style(A::AbstractMatrix{T})::Matrix{T} where {T<:Union{Int64,Float64}}
    dim = size(A)[1]
    return reshape(A, dim, dim)
end

function matrix_style(A::T)::Matrix{T} where {T<:Union{Int64,Float64}}
    return reshape([A], 1, 1)
end

"""
    get_fourth_order(CVoigt, dof)

Constructs a symmetric fourth-order tensor from a Voigt notation vector. It uses Tensors.jl package.

This function takes a Voigt notation vector `CVoigt` and the degree of freedom `dof`
to create a symmetric fourth-order tensor. The `CVoigt` vector contains components
that represent the tensor in Voigt notation, and `dof` specifies the dimension
of the tensor.

# Arguments
- `CVoigt::Matrix{Float64}`: A vector containing components of the tensor in Voigt notation.
- `dof::Int64`: The dimension of the resulting symmetric fourth-order tensor.

# Returns
- `SymmetricFourthOrderTensor{dof}`: A symmetric fourth-order tensor of dimension `dof`.

# Example
```julia
CVoigt = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
dof = 3
result = get_fourth_order(CVoigt, dof)
"""
function get_fourth_order(CVoigt::T, ::Val{DOF}) where {DOF,T}
    return fromvoigt(SymmetricTensor{4,DOF}, CVoigt)
end

"""
    find_inverse_bond_id(nlist::Vector{Vector{Int64}})

Finds the inverse of the bond id in the nlist.

# Arguments
- `nlist::Vector{Vector{Int64}}`: The nlist to find the inverse of.
# Returns
- `inverse_nlist::Vector{Dict{Int64,Int64}}`: The inverse nlist.
"""
function find_inverse_bond_id(nlist::Vector{Vector{Int64}})
    inverse_nlist = [Dict{Int64,Int64}() for _ in eachindex(nlist)]
    for iID in eachindex(nlist)
        neighbors = nlist[iID]
        for (jID, neighborID) in enumerate(neighbors)
            value = findfirst(isequal(iID), nlist[neighborID])
            # Check if neighbor ID is found in nlist, due to different horizons
            if isnothing(value)
                continue
            end

            inverse_nlist[neighborID][iID] = value
        end
    end
    return inverse_nlist
end

function get_dependent_value(datamanager::Module,
                             field_name::String,
                             parameter::Dict,
                             iID::Int64 = 1)
    dependend_value, dependent_field = is_dependent(field_name, parameter, datamanager)

    return dependend_value ?
           interpol_data(dependent_field[iID], parameter[field_name]["Data"]) :
           parameter[field_name]
end

function is_dependent(field_name::String, damage_parameter::Dict, datamanager::Module)
    if haskey(damage_parameter, field_name) && damage_parameter[field_name] isa Dict
        if !datamanager.has_key(damage_parameter[field_name]["Field"] * "NP1")
            @error "$(damage_parameter[field_name]["Field"]) does not exist for value interpolation."
            return nothing
        end
        field = datamanager.get_field(damage_parameter[field_name]["Field"], "NP1")
        return true, field
    end
    return false, nothing
end

function interpolation(x::Union{Vector{Float64},Vector{Int64}},
                       y::Union{Vector{Float64},Vector{Int64}})
    k = 3
    if length(x) <= k
        k = length(x) - 1
    end
    return Dict("spl" => Spline1D(x, y, k = k), "min" => minimum(x), "max" => maximum(x))
end

function interpol_data(x::Union{Vector{Float64},Vector{Int64},Float64,Int64},
                       values::Dict{String,Any},
                       warning_flag::Bool = true)
    if warning_flag
        if values["min"] > minimum(x)
            @warn "Interpolation value is below interpolation range. Using minimum value of dataset."
        end
        if values["max"] < maximum(x)
            @warn "Interpolation value is above interpolation range. Using maximum value of dataset."
        end
        warning_flag = false
    end
    return evaluate(values["spl"], x)
end

"""
    invert(A::Union{Matrix{Float64},Matrix{Int64}}, error_message::String="Matrix is singular")

Invert a n x n matrix. Throws an error if A is singular.

# Arguments
- A::Union{Matrix{Float64},Matrix{Int64}}: A n x n matrix.
- error_message::String="Matrix is singular": The error message returned if A is singular.

# Returns
- inverted matrix or nothing if not inverable.
"""
function invert(A::AbstractMatrix{Float64},
                error_message::String = "Matrix is singular")
    try
        return inv(smat(A))
    catch e
        @error error_message
    end
end

function determinant(A::AbstractMatrix{Float64})
    return det(smat(A))
end

function smat(A::AbstractMatrix{Float64})
    size_A = size(A)

    if size_A == (2, 2)
        return SMatrix{2,2,Float64}(A)
    elseif size_A == (3, 3)
        return SMatrix{3,3,Float64}(A)
    else
        # Fallback for other sizes, returning a copy to ensure immutability
        # or you can throw an error if unsupported sizes are not allowed
        return SMatrix{size_A[1],size_A[2],Float64}(A)
    end
end

function find_local_neighbors(nID::Int64,
                              coordinates::Union{SubArray,Matrix{Float64},Matrix{Int64}},
                              nlist::Union{Vector{Int64},SubArray{Int64}},
                              bond_horizon::Union{Float64,Int64})
    # excludes right now iID node in the coordinates list. Because it is an abritrary sublist it should be fine.
    # saving faster than recalculation?
    nlist_without_neighbor = view(nlist[nlist .!= nID], :)
    data = transpose(coordinates[nlist_without_neighbor, :])
    nnodes = length(nlist_without_neighbor)
    nhs = GridNeighborhoodSearch{size(coordinates)[2]}(search_radius = bond_horizon,
                                                       n_points = nnodes)
    initialize_grid!(nhs, data)
    neighborList = []
    foreach_neighbor(coordinates[nID, :], data, nhs, 1) do i, j, _, L
        push!(neighborList, nlist_without_neighbor[j])
    end
    return neighborList
end

"""
    progress_bar(rank::Int64, nsteps::Int64, silent::Bool)

Create a progress bar if the rank is 0. The progress bar ranges from 1 to nsteps + 1.

# Arguments
- rank::Int64: An integer to determine if the progress bar should be created.
- nsteps::Int64: The total number of steps in the progress bar.
- silent::Bool: de/activates the progress bar
# Returns
- ProgressBar or UnitRange: If rank is 0, a ProgressBar object is returned. Otherwise, a range from 1 to nsteps + 1 is returned.
"""
function progress_bar(rank::Int64, nsteps::Int64, silent::Bool)
    # Check if rank is equal to 0.
    if rank == 0 && !silent
        # If rank is 0, create and return a ProgressBar from 1 to nsteps + 1.
        return ProgressBar(1:(nsteps + 1))
    end
    # If rank is not 0, return a range from 1 to nsteps + 1.
    return 1:(nsteps + 1)
end

"""
    rotate(nodes::AbstractVector{Int64}, dof::Int64, matrix::Union{SubArray,Array{Float64,3}}, angles::SubArray, back::Bool)

Rotates the matrix.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `matrix::Union{SubArray,Array{Float64,3}}`: Matrix.
- `rot::SubArray`: Rotation tensor.
- `back::Bool`: Back.
# Returns
- `matrix::SubArray`: Matrix.
"""
function rotate(nodes::AbstractVector{Int64},
                matrix::Union{SubArray,Array{Float64,3}},
                rot::Union{SubArray,Array{Float64,3}},
                back::Bool)
    for iID in nodes
        matrix[iID, :,
               :] = rotate_second_order_tensor(rot[iID, :, :], matrix[iID, :, :],
                                               back)
    end
    return matrix
end

"""
    rotate_second_order_tensor(angles::Union{Vector{Float64},Vector{Int64}}, tensor::Matrix{Float64}, dof::Int64, back::Bool)

Rotates the second order tensor.

# Arguments
- `angles::Union{Vector{Float64},Vector{Int64}}`: Angles.
- `tensor::Matrix{Float64}`: Second order tensor.
- `dof::Int64`: Degree of freedom.
- `back::Bool`: Back.
# Returns
- `tensor::Matrix{Float64}`: Second order tensor.
"""
function rotate_second_order_tensor(R::Matrix{Float64}, tensor::Matrix{Float64}, back::Bool)
    if back
        rotation(R', tensor)
    else
        rotation(R, tensor)
    end
    return tensor
end

function rotation(R::Union{Adjoint{Float64,Matrix{Float64}},Matrix{Float64}},
                  tensor::Matrix{Float64})
    @inbounds @fastmath for m in axes(tensor, 1)
        @inbounds @fastmath for n in axes(tensor, 1)
            tmn = zero(eltype(tensor))
            @inbounds @fastmath for i in axes(tensor, 1)
                @inbounds @fastmath for j in axes(tensor, 1)
                    tmn += tensor[i, j] * R[i, m] * R[j, n]
                end
            end
            tensor[m, n] = tmn
        end
    end
end

end
