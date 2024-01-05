# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using ProgressBars
using Tensors

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
function find_active(active::Vector{Bool})
    return [i for (i, is_active) in enumerate(active) if is_active]
end

"""
    find_updatable(active::SubArray, update_list::SubArray)

Returns the indices of `active` that are true.

# Arguments
- `active::SubArray`: The vector to search in.
- `update_list::SubArray`: The vector to search in.
# Returns
- `indices::Vector`: The indices of `active` that are true.
"""
function find_updatable(active::Vector{Int64}, update_list::SubArray)
    return [active[i] for i in eachindex(active) if update_list[i] == 1]
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
    find_files_with_ending(folder_path::AbstractString, file_ending::AbstractString)

Returns a list of files in `folder_path` that end with `file_ending`.

# Arguments
- `folder_path::AbstractString`: The path to the folder.
- `file_ending::AbstractString`: The ending of the files.
# Returns
- `file_list::Vector{String}`: The list of files that end with `file_ending`.
"""
function find_files_with_ending(folder_path::AbstractString, file_ending::AbstractString)
    file_list = filter(x -> isfile(joinpath(folder_path, x)) && endswith(x, file_ending), readdir(folder_path))
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
        @error msg * " is infinite."
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
function get_fourth_order(CVoigt::Matrix{Float64}, dof::Int64)
    return fromvoigt(SymmetricTensor{4,dof}, CVoigt)
end

"""
    find_inverse_bond_id(nlist::SubArray)

Finds the inverse of the bond id in the nlist.

# Arguments
- `nlist::SubArray`: The nlist to find the inverse of.
# Returns
- `inverse_nlist::Vector{Dict{Int64,Int64}}`: The inverse nlist.
"""
function find_inverse_bond_id(nlist::SubArray)
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