
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
                append!(filenames, [outputs[output]["Output Filename"]])
            end
        end
        filenames = search_for_duplicates(filenames)

        return filenames
    end
    return []
end