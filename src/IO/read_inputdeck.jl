# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using YAML: load_file, ParserError
using ..ParameterSpec: report!, strict_mode
using ..InputDeck: read_input as read_typed_input

export read_input_deck, validate_input

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
    validate_input(params; directory = "", no_strict = false) -> (deck, input)

Validates a loaded input deck against the typed input declarations. Reports
every problem at once and aborts if there is an error. Returns the
`params["PeriLab"]` dict and the typed `PeriLabInput`.
"""
function validate_input(params::Dict; directory::AbstractString = "", no_strict::Bool = false)
    if !haskey(params, "PeriLab") || !(params["PeriLab"] isa AbstractDict) ||
       length(params["PeriLab"]) < 2
        @abort "Yaml file is not valid."
        return
    end
    deck = params["PeriLab"]
    input, ctx = read_typed_input(deck, directory;
                                  strict = strict_mode(deck; no_strict_flag = no_strict))
    report!(ctx)
    return deck, input
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
