# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Read_Input_Deck
using YAML
using DataFrames
include("../Support/Parameters/parameter_handling.jl")
export read_input_file

function read_input(filename::String)
    return YAML.load_file(filename)["PeriLab"]
end

function read_input_file(filename::String)
    params = Dict{String,Any}()
    if occursin("yaml", filename)
        @info "Read input file $filename"
        params = read_input(filename)
    elseif occursin("xml", filename)
        @info "Read input file $filename"
        #read_xml(filename = filename)
    else
        parameter = ""
        @error "Not a supported filetype  $filename"
    end
    check_key_elements(params)
    return params
end

end