# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Read_Input_Deck
using YAML
using DataFrames
include("../Support/Parameters/parameter_handling.jl")
export read_input_file

function read_input(filename::String)
    try
        return YAML.load_file(filename)["PeriLab"]
    catch
        @error "No compatible Yaml file."
        return
    end
end

function read_input_file(filename::String)
    params = Dict{String,Any}()
    if !isfile(filename)
        @error "$filename does not exist."
        return
    end
    if occursin("yaml", filename)
        @info "Read input file $filename"
        params = read_input(filename)
    else
        @error "Not a supported filetype  $filename"
        return
    end
    return check_key_elements(params)

end

end