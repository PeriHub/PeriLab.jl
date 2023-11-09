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
function merge_exodus_files(exos)
    for exo in exos
        filename = exo.file_name
        if ".0" == filename[end-1:end]
            @info "Merge output file " * filename
            Write_Exodus_Results.merge_exodus_file(filename)
        end
    end
end

function open_result_file(result_file)
    if result_file isa Exodus.ExodusDatabase
        result_file = ExodusDatabase(result_file.file_name, "rw")
    elseif result_file["file"] isa IOStream
        result_file["file"] = open(result_file["filename"], "a")
    else
        @warn "Unknown result file type"
    end
end

function close_result_file(result_file)
    if result_file isa Exodus.ExodusDatabase
        close(result_file)
    elseif result_file["file"] isa IOStream
        close(result_file["file"])
    else
        @warn "Unknown result file type"
    end
end

function close_result_files(result_files)
    for result_file in result_files
        try
            close_result_file(result_file)
        catch
            @warn "File already closed"
        end
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
            outputs[id]["Output File Type"] = "Exodus"
            push!(result_files, Write_Exodus_Results.create_result_file(filename, nnodes, dof, max_block_id, nnsets))
        elseif ".csv" == filename[end-3:end]
            if datamanager.get_rank() == 0
                push!(result_files, Write_CSV_Results.create_result_file(filename, outputs[id]))
            end
            outputs[id]["Output File Type"] = "CSV"
        end
    end

    coords = vcat(transpose(coordinates[1:nnodes, :]))
    output_frequencies = get_output_frequency(params, nsteps)
    for id in eachindex(result_files)

        if typeof(result_files[id]) == Exodus.ExodusDatabase{Int32,Int32,Int32,Float64}
            result_files[id] = Write_Exodus_Results.init_results_in_exodus(result_files[id], outputs[id], coords, block_Id[1:nnodes], Vector{Int64}(1:max_block_id), nsets)
        end
        push!(output_frequency, Dict{String,Int64}("Counter" => 0, "Output Frequency" => output_frequencies[id], "Step" => 1))

        if outputs[id]["flush_file"]
            close_result_file(result_files[id])
        end
    end

    return result_files, outputs
end

function read_input_file(filename::String)
    return Read_Input_Deck.read_input_file(filename)
end

function write_results(result_files::Vector{Any}, time::Float64, outputs::Dict, datamanager::Module)

    for id in eachindex(result_files)
        output_type = outputs[id]["Output File Type"]
        # step 1 ist the zero step?!
        output_frequency[id]["Counter"] += 1
        if output_frequency[id]["Counter"] == output_frequency[id]["Output Frequency"]
            output_frequency[id]["Step"] += 1
            nodal_outputs = Dict(key => value for (key, value) in outputs[id]["Fields"] if (!value["global_var"]))
            global_outputs = Dict(key => value for (key, value) in outputs[id]["Fields"] if (value["global_var"]))
            if outputs[id]["flush_file"]
                open_result_file(result_files[id])
            end
            if output_type == "Exodus" && length(nodal_outputs) > 0 && result_files[id] isa Exodus.ExodusDatabase
                result_files[id] = Write_Exodus_Results.write_step_and_time(result_files[id], output_frequency[id]["Step"], time)
                result_files[id] = Write_Exodus_Results.write_nodal_results_in_exodus(result_files[id], output_frequency[id]["Step"], nodal_outputs, datamanager)
            end
            if length(global_outputs) > 0
                global_values = get_global_values(global_outputs, datamanager)
                if output_type == "Exodus"
                    result_files[id] = Write_Exodus_Results.write_global_results_in_exodus(result_files[id], output_frequency[id]["Step"], global_outputs, global_values)
                end
                if output_type == "CSV"
                    Write_Exodus_Results.write_global_results_in_csv(result_files[id], global_values)
                end
            end

            if outputs[id]["flush_file"]
                close_result_file(result_files[id])
            end
            output_frequency[id]["Counter"] = 0
        end
    end
    return result_files
end

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
            global_value = calculate_block(datamanager, fieldname, calculation_type, block_id)
        elseif compute_class == "Nodeset_Data"
            node_set = get_node_set(output[varname]["compute_params"])
            global_value = calculate_nodelist(datamanager, fieldname, calculation_type, node_set)
        end

        append!(global_values, global_value)
    end
    return global_values
end

function show_block_summary(solver_options::Dict, params::Dict, datamanager::Module)
    headers = ["Block", "Material", "Damage", "Thermal", "Additive", "Density", "Horizon", "Number of Nodes"]
    df = DataFrame([header => [] for header in headers])
    # tbd
    #types = [Int64, String, String, String, String, Float64, Float64, Int64]
    #df = DataFrame([header => Vector{t}() for (header, t) in zip(headers, types)])
    #---
    block_Id = datamanager.get_field("Block_Id")
    block_list = datamanager.get_block_list()
    block_list = ["block_" * string(block) for block in block_list]

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
        push!(row, string(length(findall(x -> x == id, block_Id))))
        push!(df, row)
    end

    @info df

end

end