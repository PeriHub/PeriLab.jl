function find_files_with_ending(folder_path::AbstractString, file_ending::AbstractString)
    file_list = filter(x -> isfile(joinpath(folder_path, x)) && endswith(x, file_ending), readdir(folder_path))
    return file_list
end

function check_inf_or_nan(array, msg)
    if !isfinite(sum(array))
        @error msg * " are infinite."
        return true
    end
    return false
end

function matrix_style(A)
    """
    include a scalar or an array and reshape it to style needed for LinearAlgebra package
    """
    if length(size(A)) == 0
        A = [A]
    end
    dim = size(A)[1]
    return reshape(A, dim, dim)
end