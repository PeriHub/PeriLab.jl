# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using ProgressBars
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