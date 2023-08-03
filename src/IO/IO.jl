module IO

import .Read_Input_Deck
import .Read_Mesh
import .Write_Exodus_Results


include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")
export initialize_data

function initialize_data(filename, datamanager, comm)
    return Read_Mesh.init_data(read_input_file(filename), datamanager, comm)
end


function read_input_file(filename)
    return Read_Input_Deck.read_input_file(filename)
end

function init_write_results(params, datamanager)
    filenames = get_output_filenames(params)
    dof = datamanager.get_dof()
    coordinates = datamanager.get_field("Coordinates")'
    block_Id = datamanager.get_field("Block_Id")
    volume = datamanager.get_field("Volume")


    if dof == 2
        coords = [coordinates[:, 1], coordinates[:, 2]]
    else
        coords = [coordinates[:, 1], coordinates[:, 2], coordinates[:, 3]]
    end

    for file in filenames
        init_write_results_in_exodus(file, coords, block_Id, volume)
    end
    # coords = [
    #     1.0 0.5 0.5 1.0 0.0 0.0 0.5 1.0 0.0
    #     1.0 1.0 0.5 0.5 1.0 0.5 0.0 0.0 0.0
    #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    # ]
    return filenames
end

function write_results(params, datamanager)



end

#export read_mesh
#export load_mesh_and_distribute
end