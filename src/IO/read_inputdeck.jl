# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using YAML: load_file
using ...Parameter_Handling: validate_yaml

export read_input_file

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
        # if isa(e, YAML.ParserError)
        #     @error "Yaml Parser Error. Make sure the yaml file is valid."
        # end
        @error "Failed to read $filename."
        return nothing
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
function read_input_file(filename::String)
    params = Dict{String,Any}()
    if !isfile(filename)
        @error "$filename does not exist."
        return nothing
    end
    if !occursin("yaml", filename)
        @error "Not a supported filetype $filename"
        return nothing
    end
    @info "Read input file $filename"
    return validate_yaml(read_input(filename))
end
