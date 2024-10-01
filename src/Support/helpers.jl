# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Helpers
using PointNeighbors: GridNeighborhoodSearch, initialize_grid!, foreach_neighbor
using Tensors
using Dierckx
using ProgressBars
using LinearAlgebra
using TensorOperations
using StaticArrays
using LoopVectorization
export qdim
export check_inf_or_nan
export find_active
export get_active_update_nodes
export find_indices
export find_inverse_bond_id
export find_files_with_ending
export matrix_style
export get_fourth_order
export interpolation
export interpol_data
export progress_bar
export invert
export rotate
export fastdot


function fastdot(a, b)
    s = 0.0
    @inbounds @simd for i ∈ eachindex(a, b)
        s += a[i] * b[i]
    end
    s
end


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
    return sum(binomial(i + dof - 1, dof - 1) for i = 1:order)
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
function find_indices(vector, what)
    return findall(item -> item == what, vector)
end

"""
    find_active(active::Vector{Bool})

Returns the indices of `active` that are true.

# Arguments
- `active::Vector{Bool}`: The vector to search in.
# Returns
- `indices::Vector`: The indices of `active` that are true.
"""
function find_active(active::Union{Vector{Bool},BitVector})
    temp = zeros(Int64, sum(active))
    count::Int64 = 1
    for (i, is_active) in enumerate(active)
        if is_active
            temp[count] = i
            count += 1
        end

    end
    return temp
    #return [i for (i, is_active) in enumerate(active) if is_active]
end

"""
    get_active_update_nodes(active::Vector{Bool}, update_list::Vector{Bool}, block_nodes::Dict{Int64,Vector{Int64}}, block::Int64)

Returns the active nodes and the update nodes.

# Arguments
- `active::Vector{Bool}`: The active vector.
- `update_list::Vector{Bool}`: The update vector.
- `block_nodes::Dict{Int64,Vector{Int64}}`: The vector to search in.
- `block::Int64`: The block_id.
# Returns
- `active_nodes::Vector{Int64}`: The nodes of `active` that are true.
- `update_nodes::Vector{Int64}`: The nodes of `update` that are true.
"""
function get_active_update_nodes(
    active::Vector{Bool},
    update_list::Vector{Bool},
    block_nodes::Dict{Int64,Vector{Int64}},
    block::Int64,
)
    active_nodes::Vector{Int64} = []
    update_nodes::Vector{Int64} = []
    active_index = find_active(active[block_nodes[block]])

    if isempty(active_index)
        return active_nodes, update_nodes
    end
    active_nodes = view(block_nodes[block], active_index)
    update_index = find_active(update_list[active_nodes])
    if !isempty(update_index)
        @views update_nodes = active_nodes[update_index]
    end
    return active_nodes, update_nodes
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
    file_list = filter(
        x -> isfile(joinpath(folder_path, x)) && endswith(x, file_ending),
        readdir(folder_path),
    )
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
function check_inf_or_nan(array, msg)
    if !isfinite(sum(array))
        @error "Field ''$msg'' is infinite."
        return true
    end
    if any(isnan.(array))
        @error "Field ''$msg'' has NaN elements."
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
function matrix_style(A)

    if length(size(A)) == 0
        A = [A]
    end
    dim = size(A)[1]
    return reshape(A, dim, dim)
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
function get_fourth_order(CVoigt, dof::Int64)
    return fromvoigt(SymmetricTensor{4,dof}, CVoigt)
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

function interpolation(
    x::Union{Vector{Float64},Vector{Int64}},
    y::Union{Vector{Float64},Vector{Int64}},
)
    k = 3
    if length(x) <= k
        k = length(x) - 1
    end
    return Dict("spl" => Spline1D(x, y, k = k), "min" => minimum(x), "max" => maximum(x))
end

function interpol_data(
    x::Union{Vector{Float64},Vector{Int64},Float64,Int64},
    values::Dict{String,Any},
)
    if values["min"] > minimum(x)
        @warn "Interpolation value is below interpolation range. Using minimum value of dataset."
    end
    if values["max"] < maximum(x)
        @warn "Interpolation value is above interpolation range. Using maximum value of dataset."
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
function invert(
    A::Union{Matrix{Float64},Matrix{Int64},SubArray{Float64},SubArray{Int64},MMatrix},
    error_message::String = "Matrix is singular",
)
    try
        return inv(A[:, :])
    catch
        @error error_message
        return nothing
    end

end


function find_local_neighbors(
    nID::Int64,
    coordinates::Union{SubArray,Matrix{Float64},Matrix{Int64}},
    nlist::Union{Vector{Int64},SubArray{Int64}},
    bond_horizon::Union{Float64,Int64},
)
    # excludes right now iID node in the coordinates list. Because it is an abritrary sublist it should be fine.
    # saving faster than recalculation?
    nlist_without_neighbor = view(nlist[nlist.!=nID], :)
    data = transpose(coordinates[nlist_without_neighbor, :])
    nnodes = length(nlist_without_neighbor)
    nhs = GridNeighborhoodSearch{size(coordinates)[2]}(
        search_radius = bond_horizon,
        n_points = nnodes,
    )
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
        return ProgressBar(1:nsteps+1)
    end
    # If rank is not 0, return a range from 1 to nsteps + 1.
    return 1:nsteps+1
end


"""
    rotate(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, matrix::Union{SubArray,Array{Float64,3}}, angles::SubArray, back::Bool)

Rotates the matrix.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `matrix::Union{SubArray,Array{Float64,3}}`: Matrix.
- `rot::SubArray`: Rotation tensor.
- `back::Bool`: Back.
# Returns
- `matrix::SubArray`: Matrix.
"""
function rotate(
    nodes::Union{SubArray,Vector{Int64}},
    matrix::Union{SubArray,Array{Float64,3}},
    rot::Union{SubArray,Array{Float64,3}},
    back::Bool,
)
    for iID in nodes
        matrix[iID, :, :] =
            rotate_second_order_tensor(rot[iID, :, :], matrix[iID, :, :], back)
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
        @tensor begin
            tensor[m, n] = tensor[i, j] * R[m, i] * R[n, j]
        end
    else
        @tensor begin
            tensor[m, n] = tensor[i, j] * R[i, m] * R[j, n]
        end
    end
end

end
