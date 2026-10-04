# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export get_flush_file
export get_write_after_damage
export get_start_time
export get_end_time
export get_outputs
export get_output_frequencies
export get_output_filenames
export get_output_type
export check_for_duplicates
export output_filenames, output_frequencies, output_fieldnames

"""
    check_for_duplicates(filenames)

Check for duplicate filenames.

# Arguments
- `filenames::Vector{String}`: The filenames
"""
function check_for_duplicates(filenames::Vector{String})
    returnfilenames = []
    checked_filenames = []
    for filename in filenames
        if (filename in checked_filenames) == false
            num_same_filenames = length(findall(x -> x == filename, filenames))
            if num_same_filenames > 1
                @abort "Filename $filename is used $num_same_filenames times"
                return
            end
        end
    end
    return false
end

"""
    get_output_filenames(params::Dict, output_dir::String)

Gets the output filenames.

# Arguments
- `params::Dict`: The parameters
- `output_dir::String`: The file directory
# Returns
- `filenames::Vector{String}`: The filenames
"""
function get_output_filenames(params::Dict, output_dir::String)
    if haskey(params::Dict, "Outputs")
        filenames::Vector{String} = []
        outputs = params["Outputs"]
        for output in keys(outputs)
            output_type = get_output_type(outputs, output)
            if haskey(outputs[output], "Output Filename")
                filename = outputs[output]["Output Filename"]
                if output_type == "CSV"
                    filename = filename * ".csv"
                else
                    filename = filename * ".e"
                end
                push!(filenames, joinpath(output_dir, filename))
            end
        end
        check_for_duplicates(filenames)

        return filenames
    end
    return []
end

"""
    get_output_type(outputs::Dict, output::String)

Gets the output type.

# Arguments
- `outputs::Dict`: The outputs
- `output::String`: The output
# Returns
- `output_type::String`: The output type
"""
function get_output_type(outputs::Dict, output::String)
    if haskey(outputs[output], "Output File Type")
        return outputs[output]["Output File Type"]
    else
        @warn "No Output File Type defined for $output, defaulting to Exodus"
        return "Exodus"
    end
end

"""
    get_flush_file(outputs::Dict, output::String)

Gets the flush file.

# Arguments
- `outputs::Dict`: The outputs
- `output::String`: The output
# Returns
- `flush_file::Bool`: The flush file
"""
function get_flush_file(outputs::Dict, output::String)
    get(outputs[output], "Flush File", true)
end

"""
    get_write_after_damage(outputs::Dict, output::String)

Get the write after damage.

# Arguments
- `outputs::Dict`: The outputs
- `output::String`: The output
# Returns
- `write_after_damage::Bool`: The value
"""
function get_write_after_damage(outputs::Dict, output::String)
    get(outputs[output], "Write After Damage", false)
end

"""
    get_start_time(outputs::Dict, output::String)

Get the start_time.

# Arguments
- `outputs::Dict`: The outputs
- `output::String`: The output
# Returns
- `start_time::Float64`: The value
"""
function get_start_time(outputs::Dict, output::String)
    get(outputs[output], "Start Time", 0.0)
end

"""
    get_end_time(outputs::Dict, output::String)

Get the end_time.

# Arguments
- `outputs::Dict`: The outputs
- `output::String`: The output
# Returns
- `end_time::Float64`: The value
"""
function get_end_time(outputs::Dict, output::String)
    get(outputs[output], "End Time", Inf64)
end

"""
    get_output_fieldnames(outputs::Dict, variables::Vector{String}, computes::Vector{String}, output_type::String)

Gets the output fieldnames.

# Arguments
- `outputs::Dict`: The outputs
- `variables::Vector{String}`: The variables
- `computes::Vector{String}`: The computes
- `output_type::String`: The output type
# Returns
- `output_fieldnames::Vector{String}`: The output fieldnames
"""
function get_output_fieldnames(outputs::Dict,
                               variables::Vector{String},
                               computes::Vector{String},
                               output_type::String)
    return_outputs = []
    for output in keys(outputs)
        if !isa(outputs[output], Bool)
            @abort "Output variable $output must be set to True or False"
            return
        end
        if outputs[output]
            if output_type == "CSV"
                if output in computes
                    push!(return_outputs, [output, "Constant"])
                else
                    @warn '"' * output * '"' * " is not defined as global variable"
                end
            else
                if output in variables || output in computes
                    push!(return_outputs, [output, "Constant"])
                elseif output * "NP1" in variables
                    push!(return_outputs, [output, "NP1"])
                else
                    @warn '"' * output * '"' * " is not defined as variable"
                end
            end
        end
    end
    return return_outputs
end

"""
    get_outputs(params::Dict, variables::Vector{String}, compute_names::Vector{String})

Gets the outputs.

# Arguments
- `params::Dict`: The parameters
- `variables::Vector{String}`: The variables
- `compute_names::Vector{String}`: The compute names
# Returns
- `outputs::Dict`: The outputs
"""
function get_outputs(params::Dict, variables::Vector{String}, compute_names::Vector{String})
    num = 0
    outputs = Dict()
    if haskey(params, "Outputs")
        outputs = params["Outputs"]
        for output in keys(outputs)
            output_type = get_output_type(outputs, output)
            if (haskey(outputs[output], "Output Variables")) &&
               (length(outputs[output]["Output Variables"]) > 0)
                outputs[output]["fieldnames"] = get_output_fieldnames(outputs[output]["Output Variables"],
                                                                      variables,
                                                                      compute_names,
                                                                      output_type)
            else
                @warn "No output variables are defined for " * output * "."
            end
        end
    end
    return outputs
end

"""
    get_output_frequencies(params::Dict, nsteps::Int64)

Gets the output frequencies.

# Arguments
- `params::Dict`: The parameters
- `nsteps::Int64`: The number of steps
# Returns
- `freq::Vector{Int64}`: The output frequencies
"""
function get_output_frequencies(params::Dict, nsteps::Int64,
                                step_id::Int64)
    freq = zeros(1)
    if haskey(params::Dict, "Outputs")
        outputs = params["Outputs"]
        freq = zeros(Int64, length(keys(outputs)))
        for (id, output) in enumerate(keys(outputs))
            use_frequency = true

            if haskey(outputs[output], "Number of Output Steps")
                use_frequency = false
                if haskey(outputs[output], "Output Frequency")
                    @warn "Double output step / frequency definition. First option is used. ''Output Frequency'' is ignored."
                end
            end
            value = 1
            if use_frequency
                value = outputs[output]["Output Frequency"]
            else
                value = outputs[output]["Number of Output Steps"]
            end
            if typeof(value) == String
                value = parse(Int, split(value)[step_id])
                # @abort "Output frequency or number of output steps must be an integer."
            end
            if use_frequency
                freq[id] = value
            else
                freq[id] = Int64(ceil(nsteps / value))
            end
            if freq[id] > nsteps
                freq[id] = nsteps
            end
        end
    end

    return freq
end

# Typed functions (phase 2b). The Dict getters above stay until phase 4.

"""
    output_filenames(outputs, output_dir)

Result file names in `values(outputs)` order (`.csv` for CSV outputs, `.e`
otherwise). Aborts if a name is used twice.
"""
function output_filenames(outputs::Dict{String,OutputParams}, output_dir::String)
    filenames = String[]
    for output in values(outputs)
        extension = output.output_file_type == "CSV" ? ".csv" : ".e"
        push!(filenames, joinpath(output_dir, output.output_filename * extension))
    end
    check_for_duplicates(filenames)
    return filenames
end

"""
    output_frequencies(outputs, nsteps, step_id)

Output frequency (write every n-th step) per output, in `values(outputs)`
order. `Number of Output Steps` wins over `Output Frequency`; a value given as
a string holds one entry per solver step; the result is clamped to `nsteps`.
"""
function output_frequencies(outputs::Dict{String,OutputParams}, nsteps::Int64,
                            step_id::Int64)
    frequencies = zeros(Int64, length(outputs))
    for (id, output) in enumerate(values(outputs))
        use_frequency = output.number_of_output_steps === nothing
        if !use_frequency && output.output_frequency !== nothing
            @warn "Double output step / frequency definition. First option is used. ''Output Frequency'' is ignored."
        end
        value = use_frequency ? output.output_frequency : output.number_of_output_steps
        if value isa String
            value = parse(Int, split(value)[step_id])
        end
        frequencies[id] = use_frequency ? value : Int64(ceil(nsteps / value))
        frequencies[id] = min(frequencies[id], nsteps)
    end
    return frequencies
end

"""
    output_fieldnames(variables, field_keys, compute_names, output_type)

The selected output variables as `[name, "Constant"]` or `[name, "NP1"]`.
CSV outputs only take compute classes.
"""
function output_fieldnames(variables::Dict{String,Bool}, field_keys::Vector{String},
                           compute_names::Vector{String}, output_type::String)
    fieldnames = Vector{Vector{String}}()
    for (name, selected) in variables
        selected || continue
        if output_type == "CSV"
            if name in compute_names
                push!(fieldnames, [name, "Constant"])
            else
                @warn '"' * name * '"' * " is not defined as global variable"
            end
        elseif name in field_keys || name in compute_names
            push!(fieldnames, [name, "Constant"])
        elseif name * "NP1" in field_keys
            push!(fieldnames, [name, "NP1"])
        else
            @warn '"' * name * '"' * " is not defined as variable"
        end
    end
    return fieldnames
end
