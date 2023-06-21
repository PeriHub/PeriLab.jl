using CSV
using Logging
using YAML
using DataFrames
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
    return parameter
end


function read_input(filename::String)
    return YAML.load_file(filename)["Peridigm"]
end



