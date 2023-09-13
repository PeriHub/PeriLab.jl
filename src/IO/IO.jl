module IO
include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")
include("../Support/Parameters/parameter_handling.jl")
import .Read_Input_Deck
import .Read_Mesh
import .Write_Exodus_Results
using MPI
using TimerOutputs
export close_files
export initialize_data
export init_write_results
export write_results

function close_files(exos)
    for exo in exos
        @info "Closing output file " * exo.file_name
        close(exo)
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
    outputs = get_outputs(params, datamanager.get_all_field_keys())
    mapping = Dict{Int64,Dict{String,Vector{Any}}}()
    return_outputs = Dict{Int64,Vector{String}}()

    for id in eachindex(sort(outputs))
        result_id = 0
        mapping[id] = Dict{String,Vector{Any}}()
        for fieldname in outputs[id]
            result_id += 1
            datafield = datamanager.get_field(fieldname)
            sizedatafield = size(datafield)
            if length(sizedatafield) == 0
                #if fieldname == "Forces"
                #mapping[id]["Forces"] = [fieldname, result_id, 1, typeof(datafield[1, 1])]
                # compute class must be mapped here
                @error "No field " * fieldname * " exists."
                return
            end
            if length(sizedatafield) == 1
                mapping[id][clearNP1(fieldname)] = [fieldname, result_id, 1, typeof(datafield[1, 1])]
            else
                refDof = sizedatafield[2]
                for dof in 1:refDof
                    mapping[id][clearNP1(fieldname)*Write_Exodus_Results.get_paraviewCoordinates(dof, refDof)] = [fieldname, result_id, dof, typeof(datafield[1, 1])]
                end
            end
        end
    end
    return mapping
end

function initialize_data(filename, datamanager, comm)
    datamanager.set_rank(MPI.Comm_rank(comm))
    datamanager.set_comm(comm)
    return Read_Mesh.init_data(read_input_file(filename), datamanager, comm)
end

function init_write_results(params, datamanager)
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
    nsets = datamanager.get_nsets()
    for filename in filenames
        #if datamanager.get_rank() > 0
        filename = filename * "_" * string(datamanager.get_rank())
        # end
        push!(exos, Write_Exodus_Results.create_result_file(filename, nnodes, dof, maximum(block_Id[1:nnodes]), nnsets))
    end
    coords = vcat(transpose(coordinates[1:nnodes, :]))
    println(coords)
    outputs = get_results_mapping(params, datamanager)
    for id in eachindex(exos)
        exos[id] = Write_Exodus_Results.init_results_in_exodus(exos[id], outputs[id], coords, block_Id[1:nnodes], nsets)
    end
    return exos, outputs
end

function read_input_file(filename)
    return Read_Input_Deck.read_input_file(filename)
end

function write_results(exos, step, time, outputs, datamanager)
    for id in eachindex(exos)
        # step 1 ist the zero step?!
        exos[id] = Write_Exodus_Results.write_step_and_time(exos[id], step, time)
        exos[id] = Write_Exodus_Results.write_nodal_results_in_exodus(exos[id], step, outputs[id], datamanager)
    end
    return exos
end

end