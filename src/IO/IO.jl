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
using TimerOutputs
export close_exodus_files
export close_csv_files
export initialize_data
export init_write_results
export write_results
export merge_exodus_files
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

function close_exodus_files(exos)
    for exo in exos
        @info "Closing output file " * exo.file_name
        close(exo)
    end
end

function close_csv_files(csv_files)
    for csv in csv_files
        close(csv)
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

function get_results_mapping(params, datamanager)
    compute_names = get_computes_names(params)
    outputs = get_outputs(params, datamanager.get_all_field_keys(), compute_names)
    computes = get_computes(params, datamanager.get_all_field_keys())
    output_mapping = Dict{Int64,Dict{String,Vector{Any}}}()
    compute_mapping = Dict{Int64,Dict{}}()
    return_outputs = Dict{Int64,Vector{String}}()
    # computes = Dict("External_Displacements" => Dict("Compute Class" => "Block_Data", "Calculation Type" => "Maximum", "Block" => "block_1", "Variable" => "Displacements", "Mapping" => Dict("External_Displacementsx" => Dict("result_id" => 1, "dof" => 1), "External_Displacementsy" => Dict("result_id" => 2, "dof" => 2), "External_Displacementsz" => Dict("result_id" => 3, "dof" => 3))))
    for id in eachindex(sort(outputs))
        result_id = 0
        output_mapping[id] = Dict{String,Vector{Any}}()
        compute_mapping[id] = Dict{}()
        # computes_mapping[id] = Dict{String,Vector{Any}}()
        for fieldname in outputs[id]
            result_id += 1
            global_var = false
            compute_name = ""
            #check if fieldname occursin in array computes label
            for key in keys(computes)
                if fieldname == key
                    fieldname = computes[key]["Variable"]
                    compute_name = key
                    global_var = true
                end
            end
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
            if global_var
                if length(sizedatafield) == 1
                    computes[compute_name]["Mapping"] = Dict(compute_name => Dict("result_id" => result_id, "dof" => 1))
                else
                    refDof = sizedatafield[2]
                    computes[compute_name]["Mapping"] = Dict()
                    for dof in 1:refDof
                        field_name = compute_name * Write_Exodus_Results.get_paraviewCoordinates(dof, refDof)
                        computes[compute_name]["Mapping"][field_name] = Dict("result_id" => result_id, "dof" => dof)
                    end
                end
                compute_mapping[id][compute_name] = computes[compute_name]
            else
                if length(sizedatafield) == 1
                    output_mapping[id][clearNP1(fieldname)] = [fieldname, result_id, 1, typeof(datafield[1, 1])]
                else
                    refDof = sizedatafield[2]
                    for dof in 1:refDof
                        output_mapping[id][clearNP1(fieldname)*Write_Exodus_Results.get_paraviewCoordinates(dof, refDof)] = [fieldname, result_id, dof, typeof(datafield[1, 1])]
                    end
                end
            end
        end
    end
    return output_mapping, compute_mapping
end

function initialize_data(filename::String, datamanager::Module, comm::MPI.Comm, to::TimerOutputs.TimerOutput)

    @timeit to "MPI init data" begin
        datamanager.set_rank(MPI.Comm_rank(comm))
        datamanager.set_max_rank(MPI.Comm_size(comm))
        datamanager.set_comm(comm)
    end
    return Read_Mesh.init_data(read_input_file(filename), datamanager, comm, to)

end

function init_write_results(params::Dict, datamanager::Module, nsteps::Int64)
    filenames = get_output_filenames(params)
    if length(filenames) == 0
        @warn "No futput file or output defined"
    end
    exos = []

    nnodes = datamanager.get_nnodes()
    dof = datamanager.get_dof()
    nnsets = datamanager.get_nnsets()
    coordinates = datamanager.get_field("Coordinates")
    block_Id = datamanager.get_field("Block_Id")
    max_block_id = maximum(block_Id)
    max_block_id = find_and_set_core_value_max(datamanager.get_comm(), max_block_id)
    nsets = datamanager.get_nsets()
    for filename in filenames
        if ".e" != filename[end-1:end]
            filename = filename * ".e"
        end
        if datamanager.get_max_rank() > 1
            filename = filename * "." * string(datamanager.get_max_rank()) * "." * string(datamanager.get_rank())
        end
        push!(exos, Write_Exodus_Results.create_result_file(filename, nnodes, dof, max_block_id, nnsets))
    end
    coords = vcat(transpose(coordinates[1:nnodes, :]))
    outputs, computes = get_results_mapping(params, datamanager)
    output_frequencies = get_output_frequency(params, nsteps)
    for id in eachindex(exos)

        exos[id] = Write_Exodus_Results.init_results_in_exodus(exos[id], outputs[id], computes[id], coords, block_Id[1:nnodes], Vector{Int64}(1:max_block_id), nsets)
        push!(output_frequency, Dict{String,Int64}("Counter" => 0, "Output Frequency" => output_frequencies[id], "Step" => 1))

    end
    csv_files = []
    if datamanager.get_rank() == 0 && length(computes) > 0
        csv_files = Write_CSV_Results.create_result_file(filenames, computes)
    end

    return exos, csv_files, outputs, computes
end

function read_input_file(filename::String)
    return Read_Input_Deck.read_input_file(filename)
end

function write_results(exos, csv_files, time, outputs, computes, datamanager)
    for id in eachindex(exos)
        # step 1 ist the zero step?!
        output_frequency[id]["Counter"] += 1
        if output_frequency[id]["Counter"] == output_frequency[id]["Output Frequency"]
            output_frequency[id]["Step"] += 1
            exos[id] = Write_Exodus_Results.write_step_and_time(exos[id], output_frequency[id]["Step"], time)
            exos[id] = Write_Exodus_Results.write_nodal_results_in_exodus(exos[id], output_frequency[id]["Step"], outputs[id], datamanager)
            if length(computes) > 0
                exos[id] = Write_Exodus_Results.write_global_results_in_exodus(exos[id], csv_files[id], output_frequency[id]["Step"], computes[id], datamanager)
            end
            output_frequency[id]["Counter"] = 0
        end
    end
    return exos
end

end