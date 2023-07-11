function find_files_with_ending(folder_path::AbstractString, file_ending::AbstractString)
    file_list = filter(x -> isfile(joinpath(folder_path, x)) && endswith(x, file_ending), readdir(folder_path))
    return file_list
end