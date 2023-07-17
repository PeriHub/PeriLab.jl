module IO

import .Read_Input_Deck
import .Read_Mesh
import .Write_Results

include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")


function read_input_file(filename)
    return Read_Input_Deck.read_input_file(filename)
end

function init_data(filename, datamanager, comm)
    return Read_Mesh.init_data(read_input_file(filename), datamanager, comm)
end

#export read_mesh
#export load_mesh_and_distribute
end