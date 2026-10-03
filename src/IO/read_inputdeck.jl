# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using YAML: load_file, ParserError
using ...Parameter_Handling: validate_yaml, validate_input

export read_input_file, read_input_deck

"""
    read_input(filename::String)

Reads the input deck from a yaml file

# Arguments
- `filename::String`: The name of the yaml file
# Returns
- `params::Dict{String,Any}`: The parameters read from the yaml file
"""
function read_input(filename::String)
    try
        return load_file(filename)
    catch e
        if isa(e, ParserError)
            @abort "Yaml Parser Error. Make sure the yaml file is valid."
        end
        @abort "Failed to read $filename."
    end
end

"""
    read_input_file(filename::String)

Reads the input deck from a yaml file

# Arguments
- `filename::String`: The name of the yaml file
# Returns
- `Dict{String,Any}`: The validated parameters read from the yaml file.
"""
function read_input_file(filename::String; directory::AbstractString = dirname(filename),
                         no_strict::Bool = false)
    params = Dict{String,Any}()
    if !isfile(filename)
        @abort "$(filename) can not be found. Make sure the file exist and is readable."
        return
    end
    if !occursin("yaml", filename)
        @abort "Not a supported filetype $filename"
        return
    end
    @info "Read input file $filename"
    return validate_yaml(read_input(filename); directory = directory, no_strict = no_strict)
end

"""
    read_input_deck(filename; directory = dirname(filename), no_strict = false)

Reads and validates the input deck. Returns the deck dict (for consumers not
yet switched to typed input) and the typed `PeriLabInput`.
"""
function read_input_deck(filename::String; directory::AbstractString = dirname(filename),
                         no_strict::Bool = false)
    if !isfile(filename)
        @abort "$(filename) can not be found. Make sure the file exist and is readable."
        return
    end
    if !occursin("yaml", filename)
        @abort "Not a supported filetype $filename"
        return
    end
    @info "Read input file $filename"
    return validate_input(read_input(filename); directory = directory, no_strict = no_strict)
end
