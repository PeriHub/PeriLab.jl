# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Set_modules
export include_files
export find_module_files

"""
    find_jl_files(directory::AbstractString)

Recursively find Julia files (.jl) in a directory.

This function recursively searches for Julia source files with the ".jl" extension
in the specified directory and its subdirectories. It returns a vector of file paths
for all the found .jl files.

# Arguments
- `directory::AbstractString`: The directory in which to search for .jl files.

# Returns
A vector of strings, where each string is a file path to a .jl file found in the
specified directory and its subdirectories.

# Example
```julia
jl_files = find_jl_files("/path/to/modules")
for jl_file in jl_files
    println("Found Julia file: ", jl_file)
end
"""
function find_jl_files(directory::AbstractString)
    jl_files = Vector{String}()
    if !isdir(directory)
        @error "$directory does not exists. Modules won't be loaded accurately."
        return nothing
    end

    function find_jl_recursive(current_dir::AbstractString)
        files = readdir(current_dir)
        for file in files
            file_path = joinpath(current_dir, file)
            if isfile(file_path) && endswith(file, ".jl")
                push!(jl_files, file_path)
            elseif isdir(file_path)
                find_jl_recursive(file_path)
            end
        end
    end

    find_jl_recursive(directory)

    return jl_files
end

"""
    find_module_files(directory::AbstractString, specific::String)

Search for Julia modules containing a specific function in a given directory.

This function searches for Julia modules (files with `.jl` extension) in the specified
directory and checks if they contain a specific function. It returns a list of dictionaries
where each dictionary contains the file path and the name of the module where the specific
function is found.

# Arguments
- `directory::AbstractString`: The directory to search for Julia modules.
- `specific::String`: The name of the specific function to search for.

# Returns
An array of dictionaries, where each dictionary has the following keys:
- `"File"`: The file path to the module where the specific function is found.
- `"Module Name"`: The name of the module where the specific function is found.

# Example
```julia
result = find_module_files("/path/to/modules", "my_function")
for module_info in result
    println("Function found in module: ", module_info["Module Name"])
    println("Module file path: ", module_info["File"])
end
"""
function find_module_files(directory::AbstractString, specific::String)

    files_in_folder = find_jl_files(directory)
    module_list = []
    module_name = ""
    for filename in files_in_folder

        file = open(filename, "r")
        for line in eachline(file)
            if occursin(r"\bmodule\b", line)
                module_name = split(line)[2]
            end
            if occursin("function " * specific * "()", line)
                push!(module_list, Dict("File" => filename, "Module Name" => module_name))
                break
            end
        end
        close(file)
    end
    return module_list
end

"""
    include_files(module_list::Vector{Any})

Include files specified in a list of modules.

This function iterates over a list of modules and includes the files specified
in each module's "File" key.

# Arguments
- `module_list::Vector{Any}`: A list of modules where each module is expected to
  be a dictionary-like object with a "File" key specifying the file path.

# Examples
```julia
include_files([Dict("File" => "module1.jl"), Dict("File" => "module2.jl")])
"""
function include_files(module_list::Vector{Any})
    for mod in module_list
        include(mod["File"])
    end
end

"""
    create_module_specifics(name::String, module_list::Dict{String,AbstractString}(),specifics::Dict{String,String}(), values::Tuple)

Searches for a specific function within a list of modules and calls that function if found.

This function iterates over a list of modules specified in `module_list` and looks for a module-specific function specified in the `specifics` dictionary. If the module and function are found, it calls that function with the provided `values` tuple.

# Arguments
- `name::String`: The name to match against the module names.
- `module_list::Dict{String, AbstractString}`: A dictionary of module names mapped to abstract strings.
- `specifics::Dict{String, String}`: A dictionary specifying the module-specific function to call for each module.
- `values::Tuple`: A tuple of values to be passed as arguments to the module-specific function.

# Example
```julia
module_list = Dict("Module1" => "Module1Name", "Module2" => "Module2Name")
specifics = Dict("Module1Name" => "module1_function", "Module2Name" => "module2_function")
values = (arg1, arg2)
create_module_specifics("Module1Name", module_list, specifics, values)
"""
function create_module_specifics(name::Union{String,SubString}, module_list::Vector{Any}, specifics::Dict{String,String}, values::Tuple)
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * specifics["Name"] * "()"
        if eval(Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"] * "." * specifics["Call Function"]
            function_call = eval(Meta.parse(parse_statement))
            return function_call(values...)
        end
    end
    @error "Functionality $name not found."
    return nothing
end
"""
    create_module_specifics(name::String, module_list::Dict{String,AbstractString}(),specifics::Dict{String,String}())
    # Returns: the function itself
"""
function create_module_specifics(name::Union{String,SubString}, module_list::Vector{Any}, specifics::Dict{String,String})
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * specifics["Name"] * "()"
        if eval(Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"] * "." * specifics["Call Function"]
            function_call = eval(Meta.parse(parse_statement))
            return function_call
        end
    end
    @error "Functionality $name not found."
    return nothing
end
# only module
function create_module_specifics(name::Union{String,SubString}, module_list::Vector{Any}, get_model_name::String)
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * get_model_name * "()"
        if eval(Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"]
            module_call = eval(Meta.parse(parse_statement))
            return module_call
        end
    end
    @error "Functionality $name not found."
    return nothing
end

end