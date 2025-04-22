# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module IO
using TimerOutputs
using MPI
using DataFrames
using PrettyTables
using Logging

include("../Support/Geometry.jl")

using .Geometry

include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")
include("csv_export.jl")
include("../Compute/compute_global_values.jl")
include("../Support/Parameters/parameter_handling.jl")
include("../MPI_communication/MPI_communication.jl")
using .MPI_communication: send_single_value_from_vector, synch_responder_to_controller,
                          synch_controller_to_responder,
                          synch_controller_bonds_to_responder,
                          split_vector, synch_controller_bonds_to_responder_flattened,
                          send_vector_from_root_to_core_i, send_value,
                          find_and_set_core_value_min, find_and_set_core_value_sum,
                          find_and_set_core_value_avg, gather_values, barrier
using .Parameter_Handling
using OrderedCollections: OrderedDict
export initialize_data
export init_write_results
export get_results_mapping
export init_orientations
export write_results
export merge_exodus_files
export show_block_summary
export show_mpi_summary
global_values::Vector{Float64} = []

"""
    merge_exodus_files(result_files::Vector{Any}, output_dir::String)

Merges exodus output files

# Arguments
- `result_files::Vector{Any}`: The result files
- `output_dir::String`: The file directory
"""
function merge_exodus_files(result_files::Vector{Dict}, output_dir::String)
    for result_file in result_files
        if result_file["type"] == "Exodus"
            filename = result_file["filename"]
            @info "Merge output file " * filename
            merge_exodus_file(filename)
            filename = split(basename(filename), ".")[1] * ".e"
            new_path = joinpath(output_dir, filename)
            if abspath(filename) != abspath(new_path)
                mv(filename, new_path, force = true)
                mv("epu.log", joinpath(output_dir, "epu.log"), force = true)
            end
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
# Returns
- `true`: File is closed
- `false`: File was already closed
"""
function close_result_files(result_files::Vector{Dict})
    for result_file in result_files
        try
            close_result_file(result_file)
        catch
            @warn "File already closed"
        end
    end
    return true
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
function delete_files(result_files::Vector{Dict}, output_dir::String)
    for result_file in result_files
        if result_file["type"] == "Exodus"
            # while isfile(joinpath(output_dir, "epu.log")) == false
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
    if "NP1" == name[(end - 2):end]
        return name[1:(end - 3)]
    end
    return name
end

"""
    get_results_mapping(params::Dict, path::String, datamanager::Module)

Gets the results mapping

# Arguments
- `params::Dict`: The parameters
- `path::String`: The path
- `datamanager::Module`: The datamanager
# Returns
- `output_mapping::Dict{Int64,Dict{}}`: The results mapping
"""
function get_results_mapping(params::Dict, path::String, datamanager::Module)
    compute_names = get_computes_names(params)
    outputs = get_outputs(params, datamanager.get_all_field_keys(), compute_names)
    computes = get_computes(params, datamanager.get_all_field_keys())
    output_mapping = Dict{Int64,Dict{}}()

    for (id, output) in enumerate(keys(outputs))
        output_mapping[id] = Dict{}()
        output_mapping[id]["Fields"] = Dict{}()

        fieldnames = outputs[output]["fieldnames"]
        output_mapping[id]["flush_file"] = get_flush_file(outputs, output)
        output_mapping[id]["write_after_damage"] = get_write_after_damage(outputs, output)
        output_mapping[id]["start_time"] = get_start_time(outputs, output)
        output_mapping[id]["end_time"] = get_end_time(outputs, output)
        for fieldname in fieldnames
            compute_name = ""
            compute_params = Dict{}
            global_var = false
            nodeset = []

            for key in keys(computes)
                if fieldname[1] == key
                    fieldname[1] = computes[key]["Variable"]
                    fieldname[2] = "Constant"
                    compute_name = string(key)
                    compute_params = computes[key]
                    global_var = true
                    if computes[key]["Compute Class"] == "Node_Set_Data"
                        nodeset = computes[key]["Node Set"]
                    end
                end
            end
            # end

            datafield = datamanager.get_field(fieldname[1], fieldname[2])
            sizedatafield = size(datafield)
            if length(sizedatafield) == 0
                @error "No field " * fieldname * " exists."
                return nothing
            end

            if fieldname[1] == "State Variables"
                nstatev = length(sizedatafield) == 1 ? 1 : sizedatafield[2]
                for dof in 1:nstatev
                    output_mapping[id]["Fields"]["State_Variable_" * string(dof)] = Dict("fieldname" => fieldname[1],
                                                                                         "time" => fieldname[2],
                                                                                         "global_var" => global_var,
                                                                                         "dof" => dof,
                                                                                         "type" => typeof(datafield[1,
                                                                                                                    1]))
                end
                continue
            end

            if length(sizedatafield) == 1
                if global_var
                    output_mapping[id]["Fields"][compute_name] = Dict("fieldname" => fieldname[1],
                                                                      "time" => fieldname[2],
                                                                      "global_var" => global_var,
                                                                      "dof" => 1,
                                                                      "type" => typeof(datafield[1,
                                                                                                 1]),
                                                                      "compute_params" => compute_params,
                                                                      "nodeset" => nodeset)
                else
                    output_mapping[id]["Fields"][fieldname[1]] = Dict("fieldname" => fieldname[1],
                                                                      "time" => fieldname[2],
                                                                      "global_var" => global_var,
                                                                      "dof" => 1,
                                                                      "type" => typeof(datafield[1,
                                                                                                 1]))
                end
            elseif length(sizedatafield) == 2
                i_ref_dof = sizedatafield[2]
                for dof in 1:i_ref_dof
                    if global_var
                        output_mapping[id]["Fields"][compute_name * get_paraview_coordinates(dof,
                        i_ref_dof)] = Dict("fieldname" => fieldname[1],
                                           "time" => fieldname[2],
                                           "global_var" => global_var,
                                           "dof" => dof,
                                           "type" => typeof(datafield[1, 1]),
                                           "compute_params" => compute_params,
                                           "nodeset" => nodeset)
                    else
                        output_mapping[id]["Fields"][fieldname[1] * get_paraview_coordinates(dof,
                        i_ref_dof)] = Dict("fieldname" => fieldname[1],
                                           "time" => fieldname[2],
                                           "global_var" => global_var,
                                           "dof" => dof,
                                           "type" => typeof(datafield[1, 1]))
                    end
                end
            elseif length(sizedatafield) == 3
                i_ref_dof = sizedatafield[2]
                j_ref_dof = sizedatafield[3]
                for i_dof in 1:i_ref_dof
                    for j_dof in 1:j_ref_dof
                        if global_var
                            output_mapping[id]["Fields"][compute_name * get_paraview_coordinates(i_dof,
                            i_ref_dof) * get_paraview_coordinates(j_dof, j_ref_dof)] = Dict("fieldname" => fieldname[1],
                                                                                            "time" => fieldname[2],
                                                                                            "global_var" => global_var,
                                                                                            "i_dof" => i_dof,
                                                                                            "j_dof" => j_dof,
                                                                                            "type" => typeof(datafield[1,
                                                                                                                       1,
                                                                                                                       1]),
                                                                                            "compute_params" => compute_params,
                                                                                            "nodeset" => nodeset)
                        else
                            output_mapping[id]["Fields"][fieldname[1] * get_paraview_coordinates(i_dof,
                            i_ref_dof) * get_paraview_coordinates(j_dof, j_ref_dof)] = Dict("fieldname" => fieldname[1],
                                                                                            "time" => fieldname[2],
                                                                                            "global_var" => global_var,
                                                                                            "i_dof" => i_dof,
                                                                                            "j_dof" => j_dof,
                                                                                            "type" => typeof(datafield[1,
                                                                                                                       1,
                                                                                                                       1]))
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
function initialize_data(filename::String,
                         filedirectory::String,
                         datamanager::Module,
                         comm::MPI.Comm,
                         to::TimerOutputs.TimerOutput)
    datamanager.set_directory(filedirectory)
    @timeit to "MPI init data" begin
        datamanager.set_rank(MPI.Comm_rank(comm))
        datamanager.set_max_rank(MPI.Comm_size(comm))
        datamanager.set_comm(comm)
    end
    @timeit to "init_data" datamanager, params=init_data(read_input_file(filename),
                                                         filedirectory, datamanager, comm,
                                                         to)
    return datamanager, params
end

"""
    init_orientations(datamanager::Module)

Initialize orientations.

# Arguments
- `datamanager::Module`: The datamanager
"""
function init_orientations(datamanager::Module)
    if "Angles" in datamanager.get_all_field_keys()
        datamanager.set_rotation(true)
    end
    if "Element Angles" in datamanager.get_all_field_keys()
        datamanager.set_element_rotation(true)
    end
    rotation::Bool = datamanager.get_rotation()
    element_rotation::Bool = datamanager.get_element_rotation()

    if !rotation && !element_rotation
        return datamanager
    end
    dof = datamanager.get_dof()
    nnodes = datamanager.get_nnodes()
    orientations = datamanager.create_constant_node_field("Orientations", Float64, 3)
    rotation_tensor_N, rotation_tensor_NP1 = datamanager.create_node_field("Rotation Tensor",
                                                                           Float64,
                                                                           "Matrix", dof)
    angles = datamanager.get_field("Angles")

    for iID in 1:nnodes
        rotation_tensor_N[iID, :, :] = Geometry.rotation_tensor(angles[iID, :], dof)
        rotation_tensor_NP1[iID, :, :] = rotation_tensor_N[iID, :, :]

        if dof == 2
            orientations[iID, :] = vcat(rotation_tensor_N[iID, :, 1], 0)
        elseif dof == 3
            orientations[iID, :] = rotation_tensor_N[iID, :, 1]
        end
    end
    return datamanager
end

"""
    init_write_results(params::Dict, output_dir::String, path::String, datamanager::Module, nsteps::Int64, PERILAB_VERSION::String)

Initialize write results.

# Arguments
- `params::Dict`: The parameters
- `output_dir::String`: The output directory.
- `path::String`: The path
- `datamanager::Module`: The datamanager
- `nsteps::Int64`: The number of steps
# Returns
- `result_files::Array`: The result files
- `outputs::Dict`: The outputs
"""
function init_write_results(params::Dict,
                            output_dir::String,
                            path::String,
                            datamanager::Module,
                            PERILAB_VERSION::String)
    filenames = get_output_filenames(params, output_dir)
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

    comm = datamanager.get_comm()

    block_list = datamanager.get_block_list()
    if datamanager.get_max_rank() > 1
        all_block_list = gather_values(comm, block_list)
        if datamanager.get_rank() == 0
            block_list = unique(vcat(all_block_list...))
        end
        block_list = send_value(comm, 0, block_list)
    end

    nsets = datamanager.get_nsets()
    outputs = get_results_mapping(params, path, datamanager)

    topology = nothing
    fem_block = nothing
    num_elements = 0
    elem_global_ids = nothing
    if datamanager.fem_active()
        topology = datamanager.get_field("FE Topology")
        fem_block = datamanager.get_field("FEM Block")
        num_elements = datamanager.get_num_elements()
        num_nodes_in_topo = length(unique(topology))
        # @info nnodes, num_elements, num_nodes_in_topo
        elem_global_ids = datamanager.loc_to_glob(1:(num_elements + nnodes - num_nodes_in_topo))
        # @info elem_global_ids
    end
    for name in eachindex(nsets)
        existing_nodes = intersect(global_ids, nsets[name])
        nsets[name] = datamanager.get_local_nodes(existing_nodes)
    end

    for (id, filename) in enumerate(filenames)
        rank = datamanager.get_rank()
        max_rank = datamanager.get_max_rank()
        if ".e" == filename[(end - 1):end]
            if datamanager.get_max_rank() > 1
                filename = filename *
                           "." *
                           string(max_rank) *
                           "." *
                           get_mpi_rank_string(rank, max_rank)
            end
            outputs[id]["Output File Type"] = "Exodus"
            push!(result_files,
                  create_result_file(filename,
                                     nnodes,
                                     dof,
                                     length(block_list),
                                     nnsets,
                                     num_elements,
                                     topology))
        elseif ".csv" == filename[(end - 3):end]
            if rank == 0
                push!(result_files, create_result_file(filename, outputs[id]))
            else
                push!(result_files,
                      Dict("filename" => filename, "file" => nothing, "type" => "CSV"))
            end
            outputs[id]["Output File Type"] = "CSV"
        end
    end

    coords = vcat(transpose(coordinates[1:nnodes, :]))
    for id in eachindex(result_files)
        if result_files[id]["type"] == "Exodus"
            result_files[id]["file"] = init_results_in_exodus(result_files[id]["file"],
                                                              outputs[id],
                                                              coords,
                                                              block_Id[1:nnodes],
                                                              block_list,
                                                              nsets,
                                                              global_ids,
                                                              PERILAB_VERSION,
                                                              fem_block,
                                                              topology,
                                                              elem_global_ids)
        end

        if outputs[id]["flush_file"]
            close_result_file(result_files[id])
        end
    end

    return result_files, outputs
end

"""
    set_output_frequency(params::Dict, datamanager::Module, nsteps::Int64, step_id::Int64)

Sets the output frequency.

# Arguments
- `params::Dict`: The parameters
- `datamanager::Module`: The datamanager
- `nsteps::Int64`: The number of steps
- `step_id::Int64`: The step id
"""
function set_output_frequency(params::Dict,
                              datamanager::Module,
                              nsteps::Int64,
                              step_id::Union{Nothing,Int64} = nothing)
    output_frequencies = get_output_frequency(params, nsteps)
    if isnothing(step_id) || step_id == 1
        output_frequency = []
        for id in eachindex(output_frequencies)
            push!(output_frequency,
                  Dict{String,Int64}("Counter" => 0,
                                     "Output Frequency" => output_frequencies[id],
                                     "Step" => 1))
        end
    else
        output_frequency = datamanager.get_output_frequency()
        for id in eachindex(output_frequencies)
            output_frequency[id]["Output Frequency"] = output_frequencies[id]
            output_frequency[id]["Counter"] = 0
        end
    end
    datamanager.set_output_frequency(output_frequency)
end

"""
    write_results(result_files::Vector{Any}, time::Float64, max_damage::Float64, outputs::Dict, datamanager::Module)

Write results.

# Arguments
- `result_files::Vector{Any}`: The result files
- `time::Float64`: The time
- `max_damage::Float64`: The maximum damage
- `outputs::Dict`: The outputs
- `datamanager::Module`: The datamanager
# Returns
- `result_files::Vector{Any}`: The result files
"""
function write_results(result_files::Vector{Dict},
                       time::Float64,
                       max_damage::Float64,
                       outputs::Dict,
                       datamanager::Module)
    for id in eachindex(result_files)
        output_type = outputs[id]["Output File Type"]
        # step 1 ist the zero step?!
        if outputs[id]["write_after_damage"] && max_damage == 0.0
            continue
        end
        if outputs[id]["start_time"] > time || outputs[id]["end_time"] < time
            continue
        end
        output_frequency = datamanager.get_output_frequency()
        output_frequency[id]["Counter"] += 1
        if output_frequency[id]["Counter"] == output_frequency[id]["Output Frequency"]
            output_frequency[id]["Step"] += 1
            nodal_outputs = Dict(key => value
                                 for
                                 (key, value) in outputs[id]["Fields"]
                                 if (!value["global_var"]))
            global_outputs = Dict(key => value
                                  for
                                  (key, value) in outputs[id]["Fields"]
                                  if (value["global_var"]))
            if outputs[id]["flush_file"] &&
               ((datamanager.get_rank() == 0 && output_type == "CSV") ||
                output_type == "Exodus")
                open_result_file(result_files[id])
            end
            if output_type == "Exodus" &&
               length(nodal_outputs) > 0 &&
               result_files[id]["type"] == "Exodus"
                result_files[id]["file"] = write_step_and_time(result_files[id]["file"],
                                                               output_frequency[id]["Step"],
                                                               time)
                result_files[id]["file"] = write_nodal_results_in_exodus(result_files[id]["file"],
                                                                         output_frequency[id]["Step"],
                                                                         nodal_outputs,
                                                                         datamanager)
            end
            if length(global_outputs) > 0
                global_values = get_global_values(global_outputs, datamanager)
                if output_type == "Exodus"
                    result_files[id]["file"] = write_global_results_in_exodus(result_files[id]["file"],
                                                                              output_frequency[id]["Step"],
                                                                              global_values)
                end
                if datamanager.get_rank() == 0
                    if output_type == "CSV"
                        write_global_results_in_csv(result_files[id]["file"],
                                                    time,
                                                    global_values)
                    end
                end
            end

            if outputs[id]["flush_file"] &&
               ((datamanager.get_rank() == 0 && output_type == "CSV") ||
                output_type == "Exodus")
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
    block_ids = datamanager.get_field("Block_Id")
    block_list = datamanager.get_block_list()
    for varname in keys(sort!(OrderedDict(output)))
        compute_class = output[varname]["compute_params"]["Compute Class"]
        calculation_type = output[varname]["compute_params"]["Calculation Type"]
        fieldname = output[varname]["compute_params"]["Variable"]
        global_value = 0
        dof = 1
        if haskey(output[varname], "dof")
            dof = output[varname]["dof"]
        else
            dof = [output[varname]["i_dof"], output[varname]["j_dof"]]
        end
        if compute_class == "Block_Data"
            block = output[varname]["compute_params"]["Block"]
            if block in block_list
                block_id = block_ids[findfirst(==(block), block_list)]
            else
                error("ERROR: Block name '$block' not found. Available blocks: $(join(block_list, ", "))")
            end
            global_value, nnodes = calculate_block(datamanager, fieldname, dof,
                                                   calculation_type, block_id)
        elseif compute_class == "Node_Set_Data"
            node_set = output[varname]["nodeset"]
            node_list = datamanager.get_local_nodes(node_set)
            global_value, nnodes = calculate_nodelist(datamanager, fieldname, dof,
                                                      calculation_type, node_list)
        end
        if datamanager.get_max_rank() > 1
            global_value = find_global_core_value!(global_value, calculation_type, nnodes,
                                                   datamanager)
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
function find_global_core_value!(global_value::Union{Int64,Float64},
                                 calculation_type::String,
                                 nnodes::Int64,
                                 datamanager::Module)
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
    show_block_summary(solver_options::Dict, params::Dict, log_file::String, silent::Bool, comm::MPI.Comm, datamanager::Module)

Show block summary.

# Arguments
- `solver_options::Dict`: The solver options
- `params::Dict`: The params
- `log_file::String`: The log file
- `silent::Bool`: The silent flag
- `comm::MPI.Comm`: The comm
- `datamanager::Module`: The datamanager
"""
function show_block_summary(solver_options::Dict,
                            params::Dict,
                            log_file::String,
                            silent::Bool,
                            comm::MPI.Comm,
                            datamanager::Module)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    if size == 1
        headers = [
            "Block",
            "Material",
            "Damage",
            "Thermal",
            "Corrosion",
            "Additive",
            "Density",
            "Horizon",
            "Number of Nodes"
        ]
    else
        headers = [
            "Block",
            "Material",
            "Damage",
            "Thermal",
            "Corrosion",
            "Additive",
            "Density",
            "Horizon"
        ]
    end
    df = DataFrame([header => [] for header in headers])
    # tbd
    #types = [Int64, String, String, String, String, Float64, Float64, Int64]
    #df = DataFrame([header => Vector{t}() for (header, t) in zip(headers, types)])
    #---
    block_Id = datamanager.get_field("Block_Id")
    block_list = datamanager.get_block_list()

    for id in eachindex(block_list)
        row = [block_list[id]]
        for name in headers[2:6]
            if !(name in solver_options["Models"])
                push!(row, "")
            elseif haskey(params["Blocks"][block_list[id]], name * " Model")
                push!(row, params["Blocks"][block_list[id]][name * " Model"])
            else
                push!(row, "")
            end
        end

        for name in headers[7:8]
            if haskey(params["Blocks"][block_list[id]], name)
                push!(row, string(params["Blocks"][block_list[id]][name]))
            else
                push!(row, "")
            end
        end
        if size == 1
            # get number of nodes
            num_nodes = string(length(findall(x -> x == id, block_Id)))
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
            push!(full_df, new_row)
        end
        if !silent
            pretty_table(full_df; show_subheader = false)
            stream = Logging_module.get_log_stream(2)
            if !isnothing(stream)
                pretty_table(stream, full_df; show_subheader = false)
            end
        else
            stream = Logging_module.get_log_stream(1)
            if !isnothing(stream)
                pretty_table(stream, full_df; show_subheader = false)
            end
        end
    else
        if log_file != ""
            if !silent
                pretty_table(df; show_subheader = false)
                stream = Logging_module.get_log_stream(2)
                if !isnothing(stream)
                    pretty_table(stream, df; show_subheader = false)
                end
            else
                stream = Logging_module.get_log_stream(1)
                if !isnothing(stream)
                    pretty_table(stream, df; show_subheader = false)
                end
            end
        end
    end
end

"""
    show_mpi_summary(log_file::String, silent::Bool, comm::MPI.Comm, datamanager::Module)

Show MPI summary.

# Arguments
- `log_file::String`: The log file
- `silent::Bool`: The silent flag
- `comm::MPI.Comm`: The comm
- `datamanager::Module`: The datamanager
"""
function show_mpi_summary(log_file::String,
                          silent::Bool,
                          comm::MPI.Comm,
                          datamanager::Module)
    size = MPI.Comm_size(comm)

    if size == 1
        return size
    end

    rank = MPI.Comm_rank(comm)

    block_Id = datamanager.get_field("Block_Id")
    block_list = datamanager.get_block_list()

    block_list = datamanager.get_block_list()
    if datamanager.get_max_rank() > 1
        all_block_list = gather_values(comm, block_list)
        if datamanager.get_rank() == 0
            block_list = unique(vcat(all_block_list...))
        end
        block_list = send_value(comm, 0, block_list)
    end

    nlist = datamanager.get_nlist()
    headers = ["Rank"]
    headers = vcat(headers, block_list)
    append!(headers, ["Total"])

    df = DataFrame([header => [] for header in headers])

    row = [string(rank)]
    total = 0
    for id in 1:length(block_list)
        # get number of nodes
        # num_nodes = string(length(findall(x -> x == id, block_Id)))
        if findall(x -> x == id, block_Id) == []
            push!(row, "")
            continue
        end
        num_nodes = string(sum(length.(nlist[findall(x -> x == id, block_Id)])))
        push!(row, num_nodes)
        total += parse(Int64, num_nodes)
    end
    push!(row, string(total))
    push!(df, row)
    # Gather all DataFrames to the root process (rank 0)
    all_dfs = gather_values(comm, df)

    if rank == 0
        merged_df = vcat(all_dfs...)
        if !silent
            pretty_table(merged_df; show_subheader = false)
            stream = Logging_module.get_log_stream(2)
            if !isnothing(stream)
                pretty_table(stream, merged_df; show_subheader = false)
            end
        else
            stream = Logging_module.get_log_stream(1)
            if !isnothing(stream)
                pretty_table(stream, merged_df; show_subheader = false)
            end
        end
    else
        if log_file != ""
            if !silent
                pretty_table(df; show_subheader = false)
                stream = Logging_module.get_log_stream(2)
                if !isnothing(stream)
                    pretty_table(stream, df; show_subheader = false)
                end
            else
                stream = Logging_module.get_log_stream(1)
                if !isnothing(stream)
                    pretty_table(stream, df; show_subheader = false)
                end
            end
        end
    end
end

end
