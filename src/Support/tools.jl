# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using ProgressBars
using Tensorial
function find_files_with_ending(folder_path::AbstractString, file_ending::AbstractString)
    file_list = filter(x -> isfile(joinpath(folder_path, x)) && endswith(x, file_ending), readdir(folder_path))
    return file_list
end

function check_inf_or_nan(array, msg)
    if !isfinite(sum(array))
        @error msg * " is infinite."
        return true
    end
    return false
end
"""
include a scalar or an array and reshape it to style needed for LinearAlgebra package
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

    Parameters:
    - rank::Int64: An integer to determine if the progress bar should be created.
    - nsteps::Int64: The total number of steps in the progress bar.
    - silent::Bool: de/activates the progress bar
    Returns:
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

Constructs a symmetric fourth-order tensor from a Voigt notation vector. It uses Tensorial.jl package.

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
    return fromvoigt(SymmetricFourthOrderTensor{dof}, CVoigt)
end

function find_inverse_bond_id(nlist::SubArray)
    inverse_nlist = [Dict{Int64,Int64}() for _ in 1:length(nlist)]
    for iID in eachindex(nlist)
        for (jID, neighborID) in enumerate(nlist[iID])
            value = findfirst(isequal(iID), nlist[neighborID])
            if !isnothing(value)
                inverse_nlist[neighborID][iID] = value
            end
        end
    end
    return inverse_nlist
end