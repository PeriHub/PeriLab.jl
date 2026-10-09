# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct OutputParams
    output_filename::String = req("Output Filename")
    output_file_type::String = opt("Output File Type"; default = "Exodus",
                                   allowed = ["Exodus", "CSV"])
    output_frequency::Union{Nothing,Int64,String} = opt("Output Frequency"; default = nothing,
                                                        description = "Every n-th step; per step as \"n1 n2\"")
    number_of_output_steps::Union{Nothing,Int64,String} = opt("Number of Output Steps";
                                                              default = nothing)
    output_variables::Dict{String,Bool} = req("Output Variables")
    flush_file::Bool = opt("Flush File"; default = true)
    write_after_damage::Bool = opt("Write After Damage"; default = false)
    start_time::Float64 = opt("Start Time"; default = 0.0, quantity = :time)
    end_time::Float64 = opt("End Time"; default = Inf, quantity = :time)
    bond_export::Bool = opt("Bond Export"; default = false)
    bond_blocks::Union{Nothing,Int64,String} = opt("Bond Blocks"; default = nothing,
                                                   description = "Block id or list, e.g. \"1 3\"")
end

function check!(p::OutputParams, path::String, ctx::ParseContext)
    if p.output_frequency === nothing && p.number_of_output_steps === nothing
        add_error!(ctx, path, "\"Output Frequency\" or \"Number of Output Steps\" is required")
    end
    return nothing
end
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
