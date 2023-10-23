# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

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
    if check_element(params, "Outputs")
        filenames = []
        outputs = params["Outputs"]
        for output in keys(outputs)
            output_type = get_output_type(outputs[output])
            if check_element(outputs[output], "Output Filename")
                filename = outputs[output]["Output Filename"]
                if output_type == "CSV"
                    filename = filename * ".csv"
                else
                    filename = filename * ".e"
                end
                push!(filenames, filename)
            end
        end
        filenames = search_for_duplicates(filenames)

        return filenames
    end
    return []
end

function get_output_type(output)
    output_type = "Exouds"
    if check_element(output, "Output Type")
        output_type = output["Output Type"]
    end
    return output_type
end

function get_output_fieldnames(outputs, variables, computes, output_type)
    return_outputs = String[]
    for output in keys(outputs)
        if outputs[output]
            if output_type == "CSV"
                if output in computes
                    push!(return_outputs, output)
                else
                    @warn '"' * output * '"' * " is not defined as global variable"
                end
            else
                if output in variables || output in computes
                    push!(return_outputs, output)
                elseif output * "NP1" in variables
                    push!(return_outputs, output * "NP1")
                else
                    @warn '"' * output * '"' * " is not defined as variable"
                end
            end
        end
    end
    return return_outputs
end

function get_outputs(params, variables, compute_names)
    num = 0
    if check_element(params, "Outputs")
        outputs = params["Outputs"]
        for output in keys(outputs)
            output_type = "Exouds"
            output_type = get_output_type(outputs[output])
            if check_element(outputs[output], "Output Variables")
                outputs[output]["fieldnames"] = get_output_fieldnames(outputs[output]["Output Variables"], variables, compute_names, output_type)
            else
                @warn "No output variables are defined for " * output * "."
            end

        end
    end
    return outputs
end

function get_output_frequency(params, nsteps)

    if check_element(params, "Outputs")
        outputs = params["Outputs"]
        freq = zeros(Int64, length(keys(outputs)))
        for (id, output) in enumerate(keys(outputs))
            output_options = Dict("Output Frequency" => false, "Number of Outputs" => false)
            freq[id] = 1
            for output_option in keys(output_options)
                if check_element(outputs[output], output_option)
                    if output_options[output_option]
                        @warn "Double output step / frequency definition. First option is used. $output_option is ignored."
                        continue
                    end
                    output_options[output_option] = true
                    freq[id] = outputs[output][output_option]
                    if output_options["Number of Outputs"]
                        freq[id] = Int64(ceil(nsteps / freq[id]))
                        if freq[id] < 1
                            freq[id] = 1
                        end
                    end
                    if freq[id] > nsteps
                        freq[id] = nsteps
                    end
                end
            end

        end
    end
    return freq
end