

function get_output_filenames(params)
    if check_element(params, "Output")
        filenames = []
        outputs = params["Output"]
        for output in keys(outputs)
            if check_element(outputs[output], "Output Filename")
                append!(filenames, [outputs[output]["Output Filename"]])
            end
        end
        return filenames
    end
    return []
end