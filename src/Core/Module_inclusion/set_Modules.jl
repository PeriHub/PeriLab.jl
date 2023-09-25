module Set_modules
export include_files
export include_modules
export find_module_files

function find_jl_files(directory::AbstractString)
    jl_files = Vector{String}()
    if !isdir(directory)
        @error "$directory does not exists. Modules won't be loaded accurately."
        return jl_files
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

function find_module_files(directory::AbstractString, specific)
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

function include_files(module_list)
    for mod in module_list
        #include(pwd() * mod["File"][2:end])
        include(mod["File"])
    end
end

function include_modules(module_list)
    for m in module_list
        parse_statement = "using ." * m["Module Name"]
        eval(Meta.parse(parse_statement))
    end
end

"""
    create_module_specifics(name::String, module_list::Dict{String,AbstractString}(),
                            specifics::Dict{String,String}(), values::Tuple)

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
function create_module_specifics(name::String, module_list::Vector{Any}, specifics::Dict{String,String}, values::Tuple)
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * specifics["Name"] * "()"
        if eval(Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"] * "." * specifics["Call Function"]
            function_call = eval(Meta.parse(parse_statement))
            return function_call(values...)
        end
    end
end
end