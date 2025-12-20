# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Helpers

using PointNeighbors: GridNeighborhoodSearch, initialize_grid!, foreach_neighbor
using Meshes: Ring, Point, centroid
#using Tensors
using Dierckx: Spline1D, evaluate
using ProgressBars: ProgressBar
using LinearAlgebra: Adjoint, dot, det, norm, pinv, eigvals
using StaticArrays: MMatrix, MVector, SMatrix, @SMatrix, SVector, @SVector
using LoopVectorization
using Unitful: ustrip
using CDDLib
using Polyhedra: MixedMatHRep, polyhedron, vrep, hrep

using ..Data_Manager

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
export normalize_in_place!
export fastdot
export fast_mul!
export mat_mul!
export get_mapping
export mat_mul_transpose_mat!
export get_ring
export get_hexagon
export nearest_point_id
export get_shared_horizon
export matrix_style
export matrix_to_voigt
export voigt_to_matrix
export matrix_to_vector
export vector_to_matrix

const MAPPING_2D = @SMatrix [1 1; 2 2; 2 1]
const MAPPING_3D = @SMatrix [1 1; 2 2; 3 3; 2 3; 1 3; 1 2]

@inline function get_mapping(dof::Int64)
    if dof == 2
        return MAPPING_2D
    elseif dof == 3
        return MAPPING_3D
    else
        @error "$dof is no valid mapping option."
        return Matrix{Int64}(undef, 0, 0)
    end
end
function get_shared_horizon(id)
    horizon = Data_Manager.get_field("Shared Horizon")
    return horizon[id]
end

"""
	 matrix_to_vector(mat::Matrix{T}) where {T<:Union{Float64,Int64}}

Transforming a matrix representation in a Vector{Vector} representation.

# Arguments
- `mat::Union{Matrix{Float64},Matrix{Int64}}`: Points which form the polyhedron.
# Returns
- ``: transformed data
"""
function matrix_to_vector(mat::Matrix{T}) where {T<:Union{Float64,Int64}}
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

function matrix_diff!(s3::AbstractArray{Float64,3}, nodes::AbstractVector{Int64},
                      s2::AbstractArray{Float64,3},
                      s1::AbstractArray{Float64,3})
    @views @inbounds for iID in nodes
        @inbounds @fastmath @views for m in axes(s1, 2),
                                       n in axes(s1, 3)
            s3[iID, m, n] = s2[iID, m, n] - s1[iID, m, n]
        end
    end
end

@inline function fast_mul!(stress_NP1::AbstractMatrix{Float64},
                           C::AbstractMatrix{Float64},
                           strain_increment::AbstractMatrix{Float64},
                           stress_N::AbstractMatrix{Float64},
                           mapping::SMatrix)
    @inbounds @fastmath for m in axes(mapping, 1)
        i = mapping[m, 1]
        j = mapping[m, 2]
        sNP1 = stress_N[i, j]

        for k in axes(mapping, 1)
            ki = mapping[k, 1]
            kj = mapping[k, 2]
            sNP1 += C[m, k] * strain_increment[ki, kj]
        end

        stress_NP1[i, j] = sNP1
        # Only set symmetric component if i != j
        if i != j
            stress_NP1[j, i] = sNP1
        end
    end
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
                       factor = 1) where {T<:Union{Int64,Float64}}
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
function sub_in_place!(C::Vector{T}, A::Vector{T}, B::Vector{T}) where {T<:Number}
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        C[i] = A[i] - B[i]
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
function fastdot(a::Vector{Float64}, b::Vector{Float64}, absolute = false)
    c = zero(Float64)
    @assert length(a) == length(b)
    for i in eachindex(a, b)
        if absolute
            c += abs(a[i]) * abs(b[i])
        else
            c += a[i] * b[i]
        end
    end
    c
end
function abs!(a::Vector{T}) where {T<:Union{Int64,Float64}}
    @inbounds for i in eachindex(a)
        a[i] = abs(a[i])
    end
end
function normalize_in_place!(B::Vector{T}, A::Vector{T}) where {T<:Number}
    nrm = norm(A)
    @inbounds for i in eachindex(A)
        B[i] = A[i] / nrm
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

function find_active_nodes(active::AbstractVector{Bool},
                           active_nodes::AbstractVector{Int64},
                           nodes::AbstractVector{Int64},
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

Checks if the sum of an array is finite and has only numbers. If not, an error is raised.

# Arguments
- `array`: The array to check.
- `msg`: The error message to raise.
# Returns
- `Bool`: `true` if the scalar is finite and has only numbers, `false` otherwise.
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
"""
	check_inf_or_nan(scalar, msg)

Checks if a scalar is finite and a number. If not, an error is raised.

# Arguments
- `scalar`: The scalar to check.
- `msg`: The error message to raise.
# Returns
- `Bool`: `true` if the scalar is finite, `false` otherwise.
"""
function check_inf_or_nan(scalar::T,
                          msg::String) where {T<:Union{Int64,Float64}}
    if isnan(scalar)
        @error "Scalar Value ''$msg'' is NaN."
        return true
    end
    if !isfinite(scalar)
        @error "Scalar value ''$msg'' is infinite."
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

function create_permutation(nnodes::Int64, dof::Int64)
    perm = Vector{Int}(undef, nnodes * dof)
    idx = 1
    for d in 1:dof
        for n in 1:nnodes
            old_idx = (n - 1) * dof + d  # Row-major: node, dann dof
            perm[idx] = old_idx
            idx += 1
        end
    end
    return perm
end

function create_permutation(nodes::AbstractVector{Int64}, dof::Int64)
    perm = Vector{Int}(undef, length(nodes) * dof)
    idx = 1
    for d in 1:dof
        for n in nodes
            old_idx = (n - 1) * dof + d  # Row-major: node, dann dof
            perm[idx] = old_idx
            idx += 1
        end
    end
    return perm
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
- `Array{Float64,4}`: A symmetric fourth-order tensor of dimension `dof`.

# Example
```julia
CVoigt = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
dof = 3
result = get_fourth_order(CVoigt, dof)
```
"""
const GLOBAL_TENSOR_2D = zeros(Float64, 2, 2, 2, 2)
const GLOBAL_TENSOR_3D = zeros(Float64, 3, 3, 3, 3)

function get_fourth_order(CVoigt::AbstractMatrix{Float64}, dof::Int64)
    if dof == 2
        return voigt_to_tensor4_2d!(GLOBAL_TENSOR_2D, CVoigt)
    elseif dof == 3
        return voigt_to_tensor4_3d!(GLOBAL_TENSOR_3D, CVoigt)
    else
        @error "$dof is not a valid option. Only dof 2 and 3 are supported."
        return zeros(Float64, 0, 0, 0, 0)
    end
end

function voigt_to_tensor4_2d!(T4::Array{Float64,4}, C_voigt::AbstractMatrix{Float64})
    @inbounds begin
        T4[1, 1, 1, 1] = C_voigt[1, 1]
        T4[2, 2, 2, 2] = C_voigt[2, 2]

        T4[1, 1, 2, 2] = C_voigt[1, 2]
        T4[2, 2, 1, 1] = C_voigt[2, 1]

        c33 = C_voigt[3, 3]
        T4[1, 2, 1, 2] = T4[1, 2, 2, 1] = T4[2, 1, 1, 2] = T4[2, 1, 2, 1] = c33

        c13 = C_voigt[1, 3]
        T4[1, 1, 1, 2] = T4[1, 1, 2, 1] = c13

        c23 = C_voigt[2, 3]
        T4[2, 2, 1, 2] = T4[2, 2, 2, 1] = c23

        c31 = C_voigt[3, 1]
        T4[1, 2, 1, 1] = T4[2, 1, 1, 1] = c31

        c32 = C_voigt[3, 2]
        T4[1, 2, 2, 2] = T4[2, 1, 2, 2] = c32
    end

    return T4
end
function voigt_to_tensor4_3d!(T4::Array{Float64,4}, C_voigt::AbstractMatrix{Float64})
    @inbounds begin
        c11, c22, c33 = C_voigt[1, 1], C_voigt[2, 2], C_voigt[3, 3]
        c12, c13, c23 = C_voigt[1, 2], C_voigt[1, 3], C_voigt[2, 3]
        c44, c55, c66 = C_voigt[4, 4], C_voigt[5, 5], C_voigt[6, 6]

        T4[1, 1, 1, 1] = c11
        T4[2, 2, 2, 2] = c22
        T4[3, 3, 3, 3] = c33

        T4[1, 1, 2, 2] = T4[2, 2, 1, 1] = c12
        T4[1, 1, 3, 3] = T4[3, 3, 1, 1] = c13
        T4[2, 2, 3, 3] = T4[3, 3, 2, 2] = c23

        T4[2, 3, 2, 3] = T4[2, 3, 3, 2] = T4[3, 2, 2, 3] = T4[3, 2, 3, 2] = c44
        T4[1, 3, 1, 3] = T4[1, 3, 3, 1] = T4[3, 1, 1, 3] = T4[3, 1, 3, 1] = c55
        T4[1, 2, 1, 2] = T4[1, 2, 2, 1] = T4[2, 1, 1, 2] = T4[2, 1, 2, 1] = c66

        if C_voigt[1, 4] != 0 || C_voigt[1, 5] != 0 || C_voigt[1, 6] != 0 ||
           C_voigt[2, 4] != 0 || C_voigt[2, 5] != 0 || C_voigt[2, 6] != 0 ||
           C_voigt[3, 4] != 0 || C_voigt[3, 5] != 0 || C_voigt[3, 6] != 0 ||
           C_voigt[4, 5] != 0 || C_voigt[4, 6] != 0 || C_voigt[5, 6] != 0
            c14, c15, c16 = C_voigt[1, 4], C_voigt[1, 5], C_voigt[1, 6]
            c24, c25, c26 = C_voigt[2, 4], C_voigt[2, 5], C_voigt[2, 6]
            c34, c35, c36 = C_voigt[3, 4], C_voigt[3, 5], C_voigt[3, 6]
            c45, c46, c56 = C_voigt[4, 5], C_voigt[4, 6], C_voigt[5, 6]

            T4[1, 1, 2, 3] = T4[1, 1, 3, 2] = T4[2, 3, 1, 1] = T4[3, 2, 1, 1] = c14
            T4[1, 1, 1, 3] = T4[1, 1, 3, 1] = T4[1, 3, 1, 1] = T4[3, 1, 1, 1] = c15
            T4[1, 1, 1, 2] = T4[1, 1, 2, 1] = T4[1, 2, 1, 1] = T4[2, 1, 1, 1] = c16
            T4[2, 2, 2, 3] = T4[2, 2, 3, 2] = T4[2, 3, 2, 2] = T4[3, 2, 2, 2] = c24
            T4[2, 2, 1, 3] = T4[2, 2, 3, 1] = T4[1, 3, 2, 2] = T4[3, 1, 2, 2] = c25
            T4[2, 2, 1, 2] = T4[2, 2, 2, 1] = T4[1, 2, 2, 2] = T4[2, 1, 2, 2] = c26
            T4[3, 3, 2, 3] = T4[3, 3, 3, 2] = T4[2, 3, 3, 3] = T4[3, 2, 3, 3] = c34
            T4[3, 3, 1, 3] = T4[3, 3, 3, 1] = T4[1, 3, 3, 3] = T4[3, 1, 3, 3] = c35
            T4[3, 3, 1, 2] = T4[3, 3, 2, 1] = T4[1, 2, 3, 3] = T4[2, 1, 3, 3] = c36
            T4[2, 3, 1, 3] = T4[2, 3, 3, 1] = T4[3, 2, 1, 3] = T4[3, 2, 3, 1] = c45
            T4[1, 3, 2, 3] = T4[1, 3, 3, 2] = T4[3, 1, 2, 3] = T4[3, 1, 3, 2] = c45
            T4[2, 3, 1, 2] = T4[2, 3, 2, 1] = T4[3, 2, 1, 2] = T4[3, 2, 2, 1] = c46
            T4[1, 2, 2, 3] = T4[1, 2, 3, 2] = T4[2, 1, 2, 3] = T4[2, 1, 3, 2] = c46
            T4[1, 3, 1, 2] = T4[1, 3, 2, 1] = T4[3, 1, 1, 2] = T4[3, 1, 2, 1] = c56
            T4[1, 2, 1, 3] = T4[1, 2, 3, 1] = T4[2, 1, 1, 3] = T4[2, 1, 3, 1] = c56
        end
    end

    return T4
end

"""
	find_inverse_bond_id(nlist::BondScalarState{Int64})

Finds the inverse of the bond id in the nlist.

# Arguments
- `nlist::BondScalarState{Int64}`: The nlist to find the inverse of.
# Returns
- `inverse_nlist::Vector{Dict{Int64,Int64}}`: The inverse nlist.
"""
function find_inverse_bond_id(nlist::BondScalarState{Int64})
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

function get_dependent_value(field_name::String,
                             parameter::Dict,
                             iID::Int64 = 1)
    dependend_value, dependent_field = is_dependent(field_name, parameter)

    return dependend_value ?
           interpol_data(dependent_field[iID], parameter[field_name]["Data"]) :
           parameter[field_name]
end

function is_dependent(field_name::String, damage_parameter::Dict)
    if haskey(damage_parameter, field_name) && damage_parameter[field_name] isa Dict
        if !Data_Manager.has_key(damage_parameter[field_name]["Field"] * "NP1")
            @error "$(damage_parameter[field_name]["Field"]) does not exist for value interpolation."
            return nothing
        end
        field = Data_Manager.get_field(damage_parameter[field_name]["Field"], "NP1")
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
    if abs(determinant(A)) < 1e-10
        #  @warn "$error_message \n - rank is $(rank(A))"
        return pinv(smat(A))
    end
    return inv(smat(A))
end

function determinant(A::AbstractMatrix{Float64})
    return det(smat(A))
end

function smat(A::AbstractMatrix{Float64})
    size_A = size(A)

    if size_A == (2, 2)
        return SMatrix{2,2,Float64}(A[1, 1], A[2, 1], A[1, 2], A[2, 2])
    elseif size_A == (3, 3)
        return SMatrix{3,3,Float64}(A[1, 1], A[2, 1], A[3, 1],
                                    A[1, 2], A[2, 2], A[3, 2],
                                    A[1, 3], A[2, 3], A[3, 3])
    else
        # Fallback for other sizes, returning a copy to ensure immutability
        # or you can throw an error if unsupported sizes are not allowed
        return SMatrix{size_A[1],size_A[2],Float64}(A)
    end
end

function find_local_neighbors(nID::Int64,
                              coordinates::AbstractMatrix{Float64},
                              nlist::AbstractVector{Int64},
                              bond_horizon::T) where {T<:Union{Float64,Int64}}
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
progress_bar_stable_type(::Val{true}, nsteps::Int64) = ProgressBar(1:(nsteps + 1))  # show progress
progress_bar_stable_type(::Val{false}, nsteps::Int64) = 1:(nsteps + 1)              # silent

function progress_bar(rank::Int64, nsteps::Int64, silent::Bool)
    show_progress = rank == 0 && !silent
    return progress_bar_stable_type(Val(show_progress), nsteps)
end

"""
	rotate(nodes::AbstractVector{Int64}, dof::Int64, matrix::AbstractArray{Float64,3}, angles::SubArray, back::Bool)

Rotates the matrix.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `matrix::AbstractArray{Float64,3}`: Matrix.
- `rot::SubArray`: Rotation tensor.
- `back::Bool`: Back.
# Returns
- `matrix::SubArray`: Matrix.
"""
function rotate(nodes::AbstractVector{Int64},
                matrix::AbstractArray{Float64},
                rot::AbstractArray{Float64},
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

"""
	matrix_to_voigt(matrix)

Convert a 2x2 or 3x3 matrix to Voigt notation (6x1 vector)

# Arguments
- `matrix::Matrix{Float64}`: The matrix.
# Returns
- `voigt::Vector{Float64}`: The Voigt notation.
"""

function matrix_to_voigt(matrix::AbstractMatrix{Float64})
    if length(matrix) == 4
        @inbounds return @SVector [matrix[1, 1],
            matrix[2, 2],
            0.5 * (matrix[1, 2] + matrix[2, 1])]

    elseif length(matrix) == 9
        @inbounds return @SVector [matrix[1, 1],
            matrix[2, 2],
            matrix[3, 3],
            0.5 * (matrix[2, 3] + matrix[3, 2]),
            0.5 * (matrix[1, 3] + matrix[3, 1]),
            0.5 * (matrix[1, 2] + matrix[2, 1])]
    else
        @error "Unsupported matrix size for matrix_to_voigt"
        return nothing
    end
end

"""
	voigt_to_matrix(voigt::Union{Vector{Float64},Vector{Int64}})

Convert a Voigt notation (6x1 or 3x1 vector) to a 2x2 or 3x3 matrix

# Arguments
- `voigt::Vector{Float64}`: The Voigt notation.
# Returns
- `matrix::Matrix{Float64}`: The matrix.
"""
function voigt_to_matrix(voigt::Union{MVector,SVector,Vector})
    if length(voigt) == 3
        return @SMatrix [voigt[1] voigt[3]; voigt[3] voigt[2]]
    elseif length(voigt) == 6
        return @SMatrix [voigt[1] voigt[6] voigt[5]
                         voigt[6] voigt[2] voigt[4]
                         voigt[5] voigt[4] voigt[3]]
    else
        @error "Unsupported matrix size for voigt_to_matrix"
        return nothing
    end
end

"""
	matrix_to_vector(matrix::AbstractMatrix{Float64})

Convert a 3x3 matrix to a 6x1 vector

# Arguments
- `matrix::Matrix{Float64}`: The matrix.
# Returns
- `vector::SVector{Float64}`: The Vorigt vector.
"""
function matrix_to_vector(matrix::AbstractMatrix{Float64})
    if length(matrix) == 4
        @inbounds return @SVector [
            matrix[1, 1],
            matrix[2, 2],
            0.0,
            matrix[1, 2],
            matrix[2, 1]
        ]
    elseif length(matrix) == 9
        @inbounds return @SVector [matrix[1, 1], matrix[2, 2], matrix[3, 3],
            matrix[1, 2], matrix[2, 3], matrix[3, 1],
            matrix[2, 1], matrix[3, 2], matrix[1, 3]]
    end
end

"""
	vector_to_matrix(matrix)

Convert a 6x1 vector to a 3x3 matrix

# Arguments
- `vector::Vector{Float64}`: The vector.
# Returns
- `matrix::Matrix{Float64}`: The matrix.
"""
function vector_to_matrix(vector)
    if length(vector) == 5
        return @SMatrix [vector[1] vector[3]
                         vector[4] vector[2]]
    elseif length(vector) == 9
        return @SMatrix [vector[1] vector[4] vector[9]
                         vector[7] vector[2] vector[5]
                         vector[6] vector[8] vector[3]]
    else
        @error "Unsupported vector size for vector_to_matrix"
        return nothing
    end
end

end
