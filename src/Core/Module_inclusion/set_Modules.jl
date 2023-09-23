module Module_mania
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

function find_module_files(directory::AbstractString)
    files_in_folder = find_jl_files(directory)
    module_list = []
    module_name = ""
    for filename in files_in_folder
        file = open(filename, "r")
        for line in eachline(file)
            if occursin(r"\bmodule\b", line)
                module_name = split(line)[2]
            end
            if occursin("function name()", line)
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
        include(mod["File"])
    end
end

function include_modules(module_list)
    for m in module_list
        parse_statement = "using ." * m["Module Name"]
        eval(Meta.parse(parse_statement))
    end
end

function create_module_specifics(name, module_list, specifics, values)
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * specifics["Name"] * "()"
        if eval(Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"] * "." * specifics["Call Function"] * "$values"
            return eval(Meta.parse(parse_statement))
        end
    end
end
end