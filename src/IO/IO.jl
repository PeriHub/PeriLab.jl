# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module IO
include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")
include("csv_export.jl")
include("../Compute/compute_global_values.jl")
include("../Support/Parameters/parameter_handling.jl")
include("../MPI_communication/MPI_communication.jl")
using .Read_Input_Deck
using .Read_Mesh
using .Write_Exodus_Results
using .Write_CSV_Results
using MPI
using CSV
using Exodus
using TimerOutputs
using DataFrames
export close_result_files
export initialize_data
export init_write_results
export write_results
export merge_exodus_files
export show_block_summary
output_frequency = []
global_values = []


"""
    merge_exodus_files(result_files::Vector{Any}, filedirectory::String)

Merges exodus output files

# Arguments
- `result_files::Vector{Any}`: The result files
- `filedirectory::String`: The file directory
"""
function merge_exodus_files(result_files::Vector{Dict}, filedirectory::String)

    for result_file in result_files
        if result_file["type"] == "Exodus"
            filename = result_file["filename"]
            @info "Merge output file " * filename
            Write_Exodus_Results.merge_exodus_file(filename)
            base = basename(filename)
            filename = split(base, ".")[1] * ".e"
            mv(filename, joinpath(filedirectory, filename), force=true)
            mv("epu.log", joinpath(filedirectory, "epu.log"), force=true)
        end
    end
end

"""
    open_result_file(result_file::Dict)

Opens the result file

# Arguments
- `result_file::Dict`: The result file
"""
function open_result_file(result_file::Dict)
    if result_file["type"] == "Exodus"
        result_file["file"] = ExodusDatabase(result_file["filename"], "rw")
    elseif result_file["type"] == "CSV"
        result_file["file"] = open(result_file["filename"], "a")
    end
end

"""
    close_result_file(result_file::Dict)

Closes the result file

# Arguments
- `result_file::Dict`: The result file
"""
function close_result_file(result_file::Dict)
    if !isnothing(result_file["file"])
        close(result_file["file"])
    end
end

"""
    close_result_files(result_files::Vector{Dict})

Closes the result files

# Arguments
- `result_files::Vector{Dict}`: The result files
"""
function close_result_files(result_files::Vector{Dict})
    for result_file in result_files
        try
            close_result_file(result_file)
        catch
            @warn "File already closed"
        end
    end
end

"""
    close_result_files(result_files::Vector{Dict}, outputs::Dict{Int64,Dict{}})

Closes the result files if the flush_file flag is not set

# Arguments
- `result_files::Vector{Dict}`: The result files
- `outputs::Dict{Int64,Dict{}}`: The output settings
"""
function close_result_files(result_files::Vector{Dict}, outputs::Dict{Int64,Dict{}})
    for (id, result_file) in enumerate(result_files)
        if !outputs[id]["flush_file"]
            close_result_file(result_file)
        end
    end
end

"""
    delete_files(result_files::Vector{Dict})

Deletes the result files

# Arguments
- `result_files`: The result files
"""
function delete_files(result_files::Vector{Dict}, filedirectory::String)
    for result_file in result_files
        if result_file["type"] == "Exodus"
            # while isfile(joinpath(filedirectory, "epu.log")) == false
            #     sleep(1)
            # end
            @info "Delete output file " * result_file["filename"]
            rm(result_file["filename"])
        end
    end
end

"""
    get_file_size(result_files::Vector{Dict})

Gets the file size of the result files

# Arguments
- `result_files`: The result files
# Returns
- `total_file_size`: The total file size
"""
function get_file_size(result_files::Vector{Dict})
    total_file_size = 0
    for result_file in result_files
        file_stat = stat(result_file["filename"])  # Get file information
        total_file_size += file_stat.size  # Add the file size to the total
    end
    return total_file_size
end

"""
    clearNP1(name::String)

Clears the NP1 from the name

# Arguments
- `name::String`: The name
# Returns
- `name::String`: The cleared name
"""
function clearNP1(name::String)
    if "NP1" == name[end-2:end]
        return name[1:end-3]
    end
    return name
end

"""
    get_results_mapping(params::Dict, datamanager::Module)

Gets the results mapping

# Arguments
- `params::Dict`: The parameters
- `datamanager::Module`: The datamanager
# Returns
- `output_mapping::Dict{Int64,Dict{}}`: The results mapping
"""
function get_results_mapping(params::Dict, datamanager::Module)
    compute_names = get_computes_names(params)
    outputs = get_outputs(params, datamanager.get_all_field_keys(), compute_names)
    computes = get_computes(params, datamanager.get_all_field_keys())
    output_mapping = Dict{Int64,Dict{}}()

    for (id, output) in enumerate(keys(outputs))
        result_id = 0
        output_mapping[id] = Dict{}()
        output_mapping[id]["Fields"] = Dict{}()

        fieldnames = outputs[output]["fieldnames"]
        output_type = get_output_type(outputs, output)
        flush_file = get_flush_file(outputs, output)
        output_mapping[id]["flush_file"] = flush_file
        for fieldname in fieldnames
            result_id += 1
            compute_name = ""
            compute_params = Dict{}
            global_var = false

            for key in keys(computes)
                if fieldname == key
                    fieldname = computes[key]["Variable"]
                    compute_name = string(key)
                    compute_params = computes[key]
                    global_var = true
                end
            end
            # end

            datafield = datamanager.get_field(fieldname)
            sizedatafield = size(datafield)
            if length(sizedatafield) == 0
                #if fieldname == "Forces"
                #output_mapping[id]["Forces"] = [fieldname, result_id, 1, typeof(datafield[1, 1])]
                # compute class must be mapped here
                @error "No field " * fieldname * " exists."
                return nothing
            end

            if length(sizedatafield) == 1
                if global_var
                    output_mapping[id]["Fields"][compute_name] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => 1, "type" => typeof(datafield[1, 1]), "compute_params" => compute_params)
                else
                    output_mapping[id]["Fields"][clearNP1(fieldname)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => 1, "type" => typeof(datafield[1, 1]))
                end
            elseif length(sizedatafield) == 2
                i_ref_dof = sizedatafield[2]
                for dof in 1:i_ref_dof
                    if global_var
                        output_mapping[id]["Fields"][compute_name*Write_Exodus_Results.get_paraview_coordinates(dof, i_ref_dof)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => dof, "type" => typeof(datafield[1, 1]), "compute_params" => compute_params)
                    else
                        output_mapping[id]["Fields"][clearNP1(fieldname)*Write_Exodus_Results.get_paraview_coordinates(dof, i_ref_dof)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => dof, "type" => typeof(datafield[1, 1]))
                    end
                end
            elseif length(sizedatafield) == 3
                i_ref_dof = sizedatafield[2]
                j_ref_dof = sizedatafield[3]
                for i_dof in 1:i_ref_dof
                    for j_dof in 1:j_ref_dof
                        if global_var
                            output_mapping[id]["Fields"][compute_name*Write_Exodus_Results.get_paraview_coordinates(i_dof, i_ref_dof)*Write_Exodus_Results.get_paraview_coordinates(j_dof, j_ref_dof)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "i_dof" => i_dof, "j_dof" => j_dof, "type" => typeof(datafield[1, 1, 1]), "compute_params" => compute_params)
                        else
                            output_mapping[id]["Fields"][clearNP1(fieldname)*Write_Exodus_Results.get_paraview_coordinates(i_dof, i_ref_dof)*Write_Exodus_Results.get_paraview_coordinates(j_dof, j_ref_dof)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "i_dof" => i_dof, "j_dof" => j_dof, "type" => typeof(datafield[1, 1, 1]))
                        end
                    end
                end
            end
        end
    end
    return output_mapping
end

"""
    initialize_data(filename::String, filedirectory::String, datamanager::Module, comm::MPI.Comm, to::TimerOutputs.TimerOutput)

Initialize data.

# Arguments
- `filename::String`: The name of the input file.
- `filedirectory::String`: The directory of the input file.
- `datamanager::Module`: The datamanager
- `comm::MPI.Comm`: The MPI communicator
- `to::TimerOutputs.TimerOutput`: The TimerOutput
# Returns
- `data::Dict`: The data
"""
function initialize_data(filename::String, filedirectory::String, datamanager::Module, comm::MPI.Comm, to::TimerOutputs.TimerOutput)

    @timeit to "MPI init data" begin
        datamanager.set_rank(MPI.Comm_rank(comm))
        datamanager.set_max_rank(MPI.Comm_size(comm))
        datamanager.set_comm(comm)
    end
    return Read_Mesh.init_data(read_input_file(filename), filedirectory, datamanager, comm, to)

end

"""
    init_write_results(params::Dict, filedirectory::String, datamanager::Module, nsteps::Int64)

Initialize write results.

# Arguments
- `params::Dict`: The parameters
- `filedirectory::String`: The directory of the input file.
- `datamanager::Module`: The datamanager
- `nsteps::Int64`: The number of steps
# Returns
- `result_files::Array`: The result files
- `outputs::Dict`: The outputs
"""
function init_write_results(params::Dict, filedirectory::String, datamanager::Module, nsteps::Int64)
    filenames = get_output_filenames(params, filedirectory)
    if length(filenames) == 0
        @warn "No output file or output defined"
    end
    result_files::Vector{Dict} = []

    nnodes = datamanager.get_nnodes()
    global_ids = datamanager.loc_to_glob(1:nnodes)
    dof = datamanager.get_dof()
    nnsets = datamanager.get_nnsets()
    coordinates = datamanager.get_field("Coordinates")
    block_Id = datamanager.get_field("Block_Id")
    max_block_id = maximum(block_Id)
    max_block_id = find_and_set_core_value_max(datamanager.get_comm(), max_block_id)
    nsets = datamanager.get_nsets()
    outputs = get_results_mapping(params, datamanager)
    for name in eachindex(nsets)
        existing_nodes = intersect(global_ids, nsets[name])
        nsets[name] = datamanager.get_local_nodes(existing_nodes)
    end

    for (id, filename) in enumerate(filenames)
        rank = datamanager.get_rank()
        max_rank = datamanager.get_max_rank()
        if ".e" == filename[end-1:end]
            if datamanager.get_max_rank() > 1
                filename = filename * "." * string(max_rank) * "." * get_mpi_rank_string(rank, max_rank)
            end
            outputs[id]["Output File Type"] = "Exodus"
            push!(result_files, Write_Exodus_Results.create_result_file(filename, nnodes, dof, max_block_id, nnsets))
        elseif ".csv" == filename[end-3:end]
            if rank == 0
                push!(result_files, Write_CSV_Results.create_result_file(filename, outputs[id]))
            else
                push!(result_files, Dict("filename" => filename, "file" => nothing, "type" => "CSV"))
            end
            outputs[id]["Output File Type"] = "CSV"
        end
    end

    coords = vcat(transpose(coordinates[1:nnodes, :]))
    output_frequencies = get_output_frequency(params, nsteps)
    for id in eachindex(result_files)

        if result_files[id]["type"] == "Exodus"
            result_files[id]["file"] = Write_Exodus_Results.init_results_in_exodus(result_files[id]["file"], outputs[id], coords, block_Id[1:nnodes], Vector{Int64}(1:max_block_id), nsets, global_ids)
        end
        push!(output_frequency, Dict{String,Int64}("Counter" => 0, "Output Frequency" => output_frequencies[id], "Step" => 1))

        if outputs[id]["flush_file"]
            close_result_file(result_files[id])
        end
    end

    return result_files, outputs
end

"""
    read_input_file(filename::String)

Read input file.

# Arguments
- `filename::String`: The filename
# Returns
- `data::Dict`: The data
"""
function read_input_file(filename::String)
    return Read_Input_Deck.read_input_file(filename)
end

"""
    write_results(result_files::Vector{Any}, time::Float64, outputs::Dict, datamanager::Module)

Write results.

# Arguments
- `result_files::Vector{Any}`: The result files
- `time::Float64`: The time
- `outputs::Dict`: The outputs
- `datamanager::Module`: The datamanager
# Returns
- `result_files::Vector{Any}`: The result files
"""
function write_results(result_files::Vector{Dict}, time::Float64, outputs::Dict, datamanager::Module)

    for id in eachindex(result_files)
        output_type = outputs[id]["Output File Type"]
        # step 1 ist the zero step?!
        output_frequency[id]["Counter"] += 1
        if output_frequency[id]["Counter"] == output_frequency[id]["Output Frequency"]
            output_frequency[id]["Step"] += 1
            nodal_outputs = Dict(key => value for (key, value) in outputs[id]["Fields"] if (!value["global_var"]))
            global_outputs = Dict(key => value for (key, value) in outputs[id]["Fields"] if (value["global_var"]))
            if outputs[id]["flush_file"] && ((datamanager.get_rank() == 0 && output_type == "CSV") || output_type == "Exodus")
                open_result_file(result_files[id])
            end
            if output_type == "Exodus" && length(nodal_outputs) > 0 && result_files[id]["type"] == "Exodus"
                result_files[id]["file"] = Write_Exodus_Results.write_step_and_time(result_files[id]["file"], output_frequency[id]["Step"], time)
                result_files[id]["file"] = Write_Exodus_Results.write_nodal_results_in_exodus(result_files[id]["file"], output_frequency[id]["Step"], nodal_outputs, datamanager)
            end
            if length(global_outputs) > 0
                global_values = get_global_values(global_outputs, datamanager)
                if output_type == "Exodus"
                    result_files[id]["file"] = Write_Exodus_Results.write_global_results_in_exodus(result_files[id]["file"], output_frequency[id]["Step"], global_outputs, global_values)
                end
                if datamanager.get_rank() == 0
                    if output_type == "CSV"
                        Write_Exodus_Results.write_global_results_in_csv(result_files[id]["file"], global_values)
                    end
                end
            end

            if outputs[id]["flush_file"] && ((datamanager.get_rank() == 0 && output_type == "CSV") || output_type == "Exodus")
                close_result_file(result_files[id])
            end
            output_frequency[id]["Counter"] = 0
        end
    end

    return result_files
end

"""
    get_global_values(output::Dict, datamanager::Module)

Get global values.

# Arguments
- `output::Dict`: The output
- `datamanager::Module`: The datamanager
# Returns
- `global_values::Vector`: The global values
"""
function get_global_values(output::Dict, datamanager::Module)
    global_values = []
    for varname in keys(output)
        compute_class = output[varname]["compute_params"]["Compute Class"]
        calculation_type = output[varname]["compute_params"]["Calculation Type"]
        fieldname = output[varname]["compute_params"]["Variable"]
        global_value = 0
        if compute_class == "Block_Data"
            block = output[varname]["compute_params"]["Block"]
            block_id = parse(Int, block[7:end])
            global_value, nnodes = calculate_block(datamanager, fieldname, calculation_type, block_id)
        elseif compute_class == "Nodeset_Data"
            node_set = get_node_set(output[varname]["compute_params"])
            global_value, nnodes = calculate_nodelist(datamanager, fieldname, calculation_type, node_set)
        end

        if datamanager.get_max_rank() > 1
            for iID in eachindex(global_value)
                global_value[iID] = find_global_core_value!(global_value[iID], calculation_type, nnodes, datamanager)
            end
        end
        append!(global_values, global_value)
    end
    return global_values
end

"""
    find_global_core_value!(global_value::Union{Int64,Float64}, calculation_type::String, nnodes::Int64, datamanager::Module)

Find global core value.

# Arguments
- `global_value::Union{Int64,Float64}`: The global value
- `calculation_type::String`: The calculation type
- `nnodes::Int64`: The number of nodes
- `datamanager::Module`: The datamanager
# Returns
- `global_value::Union{Int64,Float64}`: The global value
"""
function find_global_core_value!(global_value::Union{Int64,Float64}, calculation_type::String, nnodes::Int64, datamanager::Module)
    comm = datamanager.get_comm()
    if calculation_type == "Sum"
        return find_and_set_core_value_sum(comm, global_value)
    elseif calculation_type == "Maximum"
        return find_and_set_core_value_max(comm, global_value)
    elseif calculation_type == "Minimum"
        return find_and_set_core_value_min(comm, global_value)
    elseif calculation_type == "Average"
        return find_and_set_core_value_avg(comm, global_value, nnodes)
    else
        @warn "Unknown calculation type $calculation_type"
        return 0
    end
end

"""
    get_mpi_rank_string(rank::Int64, max_rank::Int64)

Get MPI rank string.

# Arguments
- `value::Int64`: The rank
- `max_rank::Int64`: The max rank
# Returns
- `result::String`: The result
"""
function get_mpi_rank_string(rank::Int64, max_rank::Int64)
    max_rank_length::Int64 = length(string(max_rank))
    rank_length::Int64 = length(string(rank))
    return "0"^(max_rank_length - rank_length) * string(rank)
end

"""
    show_block_summary(solver_options::Dict, params::Dict, datamanager::Module)

Show block summary.

# Arguments
- `solver_options::Dict`: The solver options
- `params::Dict`: The params
- `datamanager::Module`: The datamanager
"""
function show_block_summary(solver_options::Dict, params::Dict, comm::MPI.Comm, datamanager::Module)
    headers = ["Block", "Material", "Damage", "Thermal", "Additive", "Density", "Horizon", "Number of Nodes"]
    df = DataFrame([header => [] for header in headers])
    # tbd
    #types = [Int64, String, String, String, String, Float64, Float64, Int64]
    #df = DataFrame([header => Vector{t}() for (header, t) in zip(headers, types)])
    #---
    block_Id = datamanager.get_field("Block_Id")
    block_list = datamanager.get_block_list()
    block_list = ["block_" * string(block) for block in block_list]

    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    for id in eachindex(block_list)
        row = [block_list[id]]
        for name in headers[2:end-3]
            if !solver_options[name*" Models"]
                push!(row, "")
            elseif haskey(params["Blocks"][block_list[id]], name * " Model")
                push!(row, params["Blocks"][block_list[id]][name*" Model"])
            else
                push!(row, "")
            end
        end
        for name in headers[end-2:end-1]
            if haskey(params["Blocks"][block_list[id]], name)
                push!(row, string(params["Blocks"][block_list[id]][name]))
            else
                push!(row, "")
            end
        end
        # get number of nodes
        num_nodes = string(length(findall(x -> x == id, block_Id)))
        if size > 1
            push!(row, num_nodes * "($rank)")
        else
            push!(row, num_nodes)
        end
        push!(df, row)
    end

    # Gather all DataFrames to the root process (rank 0)
    all_dfs = gather_values(comm, df)

    if rank == 0 && size > 1
        merged_df = vcat(all_dfs...)
        blocks = unique(merged_df.Block)
        full_df = DataFrame()
        for block in blocks
            block_rows = filter(row -> row.Block == block, merged_df)
            new_row = block_rows[1, :]
            for row in eachrow(block_rows[2:end, :])
                if row["Number of Nodes"][1] != "0"
                    new_row["Number of Nodes"] = new_row["Number of Nodes"] * ", " * row["Number of Nodes"]
                end
            end
            push!(full_df, new_row)
        end
        @info full_df
    else
        @info df
    end

end

end