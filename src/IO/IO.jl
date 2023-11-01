# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module IO
include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")
include("csv_export.jl")
include("../Support/Parameters/parameter_handling.jl")
include("../MPI_communication/MPI_communication.jl")
using .Read_Input_Deck
using .Read_Mesh
using .Write_Exodus_Results
using .Write_CSV_Results
using MPI
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
function merge_exodus_files(exos)
    for exo in exos
        filename = exo.file_name
        if ".0" == filename[end-1:end]
            @info "Merge output file " * filename
            Write_Exodus_Results.merge_exodus_file(filename)
        end
    end
end

function close_result_files(result_files)
    for result_file in result_files
        # @info "Closing output file " * exo.file_name
        close(result_file)
    end
end

function delete_files(exos)
    for exo in exos
        @info "Delete output file " * exo.file_name
        rm(exo.file_name)
    end
end

function get_file_size(exos)
    total_file_size = 0
    for exo in exos
        file_stat = stat(exo.file_name)  # Get file information
        total_file_size += file_stat.size  # Add the file size to the total
    end
    return total_file_size
end

function clearNP1(name)
    if "NP1" == name[end-2:end]
        return name[1:end-3]
    end
    return name
end

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
        output_type = get_output_type(outputs[output])
        for fieldname in fieldnames
            result_id += 1
            compute_name = ""
            compute_params = Dict{}
            global_var = false

            for key in keys(computes)
                if fieldname == key
                    fieldname = computes[key]["Variable"]
                    compute_name = key
                    compute_params = computes[key]
                    global_var = true
                end
            end
            # end

            if datamanager.field_array_type[fieldname]["Type"] == "Matrix"
                @warn "Matrix types not supported in Exodus export"
                continue
            end
            datafield = datamanager.get_field(fieldname)
            sizedatafield = size(datafield)
            if length(sizedatafield) == 0
                #if fieldname == "Forces"
                #output_mapping[id]["Forces"] = [fieldname, result_id, 1, typeof(datafield[1, 1])]
                # compute class must be mapped here
                @error "No field " * fieldname * " exists."
                return
            end

            if length(sizedatafield) == 1
                if global_var
                    output_mapping[id]["Fields"][compute_name] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => 1, "type" => typeof(datafield[1, 1]), "compute_params" => compute_params)
                else
                    output_mapping[id]["Fields"][clearNP1(fieldname)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => 1, "type" => typeof(datafield[1, 1]))
                end
            else
                refDof = sizedatafield[2]
                for dof in 1:refDof
                    if global_var
                        output_mapping[id]["Fields"][compute_name*Write_Exodus_Results.get_paraviewCoordinates(dof, refDof)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => dof, "type" => typeof(datafield[1, 1]), "compute_params" => compute_params)
                    else
                        output_mapping[id]["Fields"][clearNP1(fieldname)*Write_Exodus_Results.get_paraviewCoordinates(dof, refDof)] = Dict("fieldname" => fieldname, "global_var" => global_var, "result_id" => result_id, "dof" => dof, "type" => typeof(datafield[1, 1]))
                    end
                end
            end
        end
    end
    return output_mapping
end

function initialize_data(filename::String, filedirectory::String, datamanager::Module, comm::MPI.Comm, to::TimerOutputs.TimerOutput)

    @timeit to "MPI init data" begin
        datamanager.set_rank(MPI.Comm_rank(comm))
        datamanager.set_max_rank(MPI.Comm_size(comm))
        datamanager.set_comm(comm)
    end
    return Read_Mesh.init_data(read_input_file(filename), filedirectory, datamanager, comm, to)

end

function init_write_results(params::Dict, filedirectory::String, datamanager::Module, nsteps::Int64)
    filenames = get_output_filenames(params, filedirectory)
    if length(filenames) == 0
        @warn "No futput file or output defined"
    end
    result_files = []

    nnodes = datamanager.get_nnodes()
    dof = datamanager.get_dof()
    nnsets = datamanager.get_nnsets()
    coordinates = datamanager.get_field("Coordinates")
    block_Id = datamanager.get_field("Block_Id")
    max_block_id = maximum(block_Id)
    max_block_id = find_and_set_core_value_max(datamanager.get_comm(), max_block_id)
    nsets = datamanager.get_nsets()
    outputs = get_results_mapping(params, datamanager)

    for (id, filename) in enumerate(filenames)

        if ".e" == filename[end-1:end]
            if datamanager.get_max_rank() > 1
                filename = filename * "." * string(datamanager.get_max_rank()) * "." * string(datamanager.get_rank())
            end
            outputs[id]["Output Type"] = "Exodus"
            push!(result_files, Write_Exodus_Results.create_result_file(filename, nnodes, dof, max_block_id, nnsets))
        elseif ".csv" == filename[end-3:end]
            if datamanager.get_rank() == 0
                push!(result_files, Write_CSV_Results.create_result_file(filename, outputs[id]))
            end
            outputs[id]["Output Type"] = "CSV"
        end
    end

    coords = vcat(transpose(coordinates[1:nnodes, :]))
    output_frequencies = get_output_frequency(params, nsteps)
    for id in eachindex(result_files)

        if typeof(result_files[id]) == Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}
            result_files[id] = Write_Exodus_Results.init_results_in_exodus(result_files[id], outputs[id], coords, block_Id[1:nnodes], Vector{Int64}(1:max_block_id), nsets)
        end
        push!(output_frequency, Dict{String,Int64}("Counter" => 0, "Output Frequency" => output_frequencies[id], "Step" => 1))

    end

    return result_files, outputs
end

function read_input_file(filename::String)
    return Read_Input_Deck.read_input_file(filename)
end

function write_results(result_files, time, outputs, datamanager)
    for id in eachindex(result_files)
        # step 1 ist the zero step?!
        output_frequency[id]["Counter"] += 1
        if output_frequency[id]["Counter"] == output_frequency[id]["Output Frequency"]
            output_frequency[id]["Step"] += 1
            nodal_outputs = Dict(key => value for (key, value) in outputs[id]["Fields"] if (!value["global_var"]))
            global_outputs = Dict(key => value for (key, value) in outputs[id]["Fields"] if (value["global_var"]))
            output_type = outputs[id]["Output Type"]
            if output_type == "Exodus" && length(nodal_outputs) > 0 && typeof(result_files[id]) == Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}
                result_files[id] = Write_Exodus_Results.write_step_and_time(result_files[id], output_frequency[id]["Step"], time)
                result_files[id] = Write_Exodus_Results.write_nodal_results_in_exodus(result_files[id], output_frequency[id]["Step"], nodal_outputs, datamanager)
            end
            if length(global_outputs) > 0
                result_files[id] = Write_Exodus_Results.write_global_results_in_exodus(result_files[id], output_frequency[id]["Step"], global_outputs, output_type, datamanager)
            end

            output_frequency[id]["Counter"] = 0
        end
    end
    return result_files
end

function show_block_summary(solver_options::Dict, params::Dict, datamanager::Module)

    headers = ["Block", "Material", "Damage", "Thermal", "Additive"]
    df = DataFrame([header => [] for header in headers])

    block_list = datamanager.get_block_list()
    block_list = ["block_" * string(block) for block in block_list]

    for id in eachindex(block_list)
        row = [block_list[id]]
        for name in headers[2:end]
            if !solver_options[name*" Models"]
                push!(row, "")
            elseif haskey(params["Blocks"][block_list[id]], name * " Model")
                push!(row, params["Blocks"][block_list[id]][name*" Model"])
            else
                push!(row, "")
            end
        end
        push!(df, row)
    end

    @info df

end

end