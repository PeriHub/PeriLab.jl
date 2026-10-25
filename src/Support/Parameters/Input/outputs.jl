# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const _EXODUS = Dict("Output File Type" => ["Exodus"])

@params "A result file with the fields and compute classes to write and when to write them" struct OutputParams
    "Name of the result file without extension"
    output_filename::String = req("Output Filename")
    "File format; CSV holds only the compute class values"
    output_file_type::String = opt("Output File Type"; default = "Exodus",
                                   allowed = ["Exodus", "CSV"])
    "Field and compute class names to write (true) or skip (false)"
    output_variables::Dict{String,Bool} = req("Output Variables")
    "Writes every n-th step; per solver step as \"n1 n2\"; Number of Output Steps takes precedence"
    output_frequency::Union{Nothing,Int64,String} = opt("Output Frequency"; default = nothing)
    "Number of outputs, evenly spread over the solver steps; per solver step as \"n1 n2\""
    number_of_output_steps::Union{Nothing,Int64,String} = opt("Number of Output Steps";
                                                              default = nothing)
    "No output before this time"
    start_time::Float64 = opt("Start Time"; default = 0.0, quantity = :time)
    "No output after this time"
    end_time::Float64 = opt("End Time"; default = Inf, quantity = :time)
    "Starts writing only once damage has occurred"
    write_after_damage::Bool = opt("Write After Damage"; default = false)
    "Closes the file after every output so it can be read while the simulation runs"
    flush_file::Bool = opt("Flush File"; default = true)
    "Writes the bonds as line elements so bond fields can be visualized; switched on automatically when a bond field is requested"
    bond_export::Bool = opt("Bond Export"; default = false, when = _EXODUS)
    "Blocks whose bonds are exported: a block id or ids separated by spaces, e.g. \"1 3\"; all blocks if not given"
    bond_blocks::Union{Nothing,Int64,
                       String} = opt("Bond Blocks"; default = nothing,
                                     when = _EXODUS)
    "Writes every n-th step of a fallback solver; per solver step as \"n1 n2\""
    fallback_output_frequency::Union{Nothing,Int64,String} = opt("Fallback Output Frequency";
                                                                 default = nothing)
    "Number of outputs over the steps of a fallback solver; per solver step as \"n1 n2\""
    fallback_number_of_output_steps::Union{Nothing,Int64,String} = opt("Fallback Number of Output Steps";
                                                                       default = nothing)
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

"A value given per solver step as a string (\"n1 n2\") is the entry of `step_id`."
step_value(value::Int64, step_id::Int64) = value
step_value(value::String, step_id::Int64) = parse(Int, split(value)[step_id])

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
        value = step_value(use_frequency ? output.output_frequency :
                           output.number_of_output_steps, step_id)
        frequencies[id] = use_frequency ? value : Int64(ceil(nsteps / value))
        frequencies[id] = min(frequencies[id], nsteps)
    end
    return frequencies
end

"""
    fallback_output_frequencies(outputs, nsteps, step_id, frequencies, counters, dt, fallback_dt)

Output frequencies and counters, in `values(outputs)` order, when a fallback
solver with `nsteps` steps of `fallback_dt` takes over from a solver with steps
of `dt` and the given `frequencies` and `counters`. `Fallback Number of Output
Steps` (which wins) and `Fallback Output Frequency` set the frequency like their
counterparts; without them an output keeps its interval in time. The counter
keeps the time since the last output. Frequencies are clamped to `nsteps`.
"""
function fallback_output_frequencies(outputs::Dict{String,OutputParams}, nsteps::Int64,
                                     step_id::Int64, frequencies::Vector{Int64},
                                     counters::Vector{Int64}, dt::Float64,
                                     fallback_dt::Float64)
    ratio = dt / fallback_dt
    new_frequencies = zeros(Int64, length(outputs))
    new_counters = zeros(Int64, length(outputs))
    for (id, output) in enumerate(values(outputs))
        frequency = if output.fallback_number_of_output_steps !== nothing
            Int64(ceil(nsteps / step_value(output.fallback_number_of_output_steps, step_id)))
        elseif output.fallback_output_frequency !== nothing
            step_value(output.fallback_output_frequency, step_id)
        else
            max(1, round(Int64, frequencies[id] * ratio))
        end
        new_frequencies[id] = min(frequency, nsteps)
        new_counters[id] = min(round(Int64, counters[id] * ratio), new_frequencies[id] - 1)
    end
    return new_frequencies, new_counters
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
