module IO

import .Read_Input_Deck
import .Read_Mesh
import .Write_Exodus_Results

include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")
include("../Support/Parameters/parameter_handling.jl")
export initialize_data
export init_write_results
function initialize_data(filename, datamanager, comm)
    return Read_Mesh.init_data(read_input_file(filename), datamanager, comm)
end

function read_input_file(filename)
    return Read_Input_Deck.read_input_file(filename)
end

function init_write_results(params, datamanager)
    filenames = get_output_filenames(params)
    exos = []
    nnodes = datamanager.get_nnodes()
    dof = datamanager.get_dof()
    nnsets = datamanager.get_nnsets()
    coordinates = datamanager.get_field("Coordinates")'
    block_Id = datamanager.get_field("Block_Id")
    nsets = datamanager.get_nsets()
    for filename in filenames
        append!(exos, Write_Exodus_Results.create_result_file(filename, nnodes, dof, maximum(block_Id), length(nsets)))
    end
    outputs = get_outputs(params, datamanager.get_all_field_keys())
    coords = vcat(transpose(coordinates))
    for num in 1:length(exos)
        exos[i] = init_results_in_exodus(exos[i], outputs[i], coords, block_Id, volume)
    end

    return exos
end

function write_results(params, datamanager)



end

#export read_mesh
#export load_mesh_and_distribute
end