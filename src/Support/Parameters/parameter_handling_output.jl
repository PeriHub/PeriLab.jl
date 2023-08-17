
function search_for_duplicates(filenames)
    returnfilenames = []
    checked_filenames = []
    for filename in filenames
        if (filename in checked_filenames) == false
            num_same_filenames = length(findall(x -> x == filename, filenames))
            if num_same_filenames > 1
                for i = 1:num_same_filenames
                    push!(returnfilenames, filename * "_" * string(i))
                end
                push!(checked_filenames, filename)
            else
                push!(returnfilenames, filename)
            end
        end
    end
    return returnfilenames
end
function get_output_filenames(params)
    if check_element(params, "Output")
        filenames = []
        outputs = params["Output"]
        for output in keys(outputs)
            if check_element(outputs[output], "Output Filename")
                push!(filenames, outputs[output]["Output Filename"])
            end
        end
        filenames = search_for_duplicates(filenames)

        return filenames
    end
    return []
end

function get_output_variables(outputs, variables)
    return_outputs = []
    for output in keys(outputs)
        if outputs[output]
            if (output in variables) || output * "NP1" in variables
                push!(return_outputs, output)
            else
                @warn '"' * output * '"' * " is not defined as variable"
            end
        end
    end
    return return_outputs
end

function get_outputs(params, variables)
    return_outputs = Dict{Int64,Any}()
    num = 0
    if check_element(params, "Output")
        outputs = params["Output"]
        for output in keys(outputs)
            if check_element(outputs[output], "Output Variables")
                num += 1
                return_outputs[num] = get_output_variables(outputs[output]["Output Variables"], variables)
            else
                @warn "No output variables are defined for " * output * "."
            end

        end
    end
    return return_outputs
end