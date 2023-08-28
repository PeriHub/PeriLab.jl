module Read_Input_Deck
using YAML
using DataFrames
include("../Support/Parameters/parameter_handling.jl")

export read_input_file

function read_input(filename::String)
    return YAML.load_file(filename)["PeriLab"]
end

function read_input_file(filename::String)
    if occursin("yaml", filename)
        println("Load  $filename")
        @info "Read file $filename"
        parameter = read_input(filename)
    elseif occursin("xml", filename)
        @info "Read file $filename"
        #read_xml(filename = filename)
    else
        parameter = ""
        @error "Not a supported filetype  $filename"
    end
    check_key_elements(parameter)
    return parameter
end



end