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
    coordinates = datamanager.get_field("Coordinates")

    block_Id = datamanager.get_field("Block_Id")
    nsets = datamanager.get_nsets()
    for filename in filenames
        push!(exos, Write_Exodus_Results.create_result_file(filename, nnodes, dof, maximum(block_Id), nnsets))
    end

    coords = vcat(transpose(coordinates))
    outputs, mapping = get_results_mapping(params, datamanager)
    for i in eachindex(exos)
        exos[i] = Write_Exodus_Results.init_results_in_exodus(exos[i], outputs[i], coords, block_Id, nsets)
    end

    return exos, mapping
end

function write_results(params, datamanager)



end

function get_results_mapping(params, datamanager)
    outputs = get_outputs(params, datamanager.get_all_field_keys())
    mapping = Dict{Int64,Dict{String,String}}()
    return_outputs = Dict{Int64,Vector{String}}()
    for id in eachindex(sort(outputs))
        return_outputs[id] = String[]
        mapping[id] = Dict{String,Vector{String}}()
        for fieldname in outputs[id]
            datafield = datamanager.get_field(fieldname)
            sizedatafield = size(datafield)
            if length(sizedatafield) == 0
                @error "No field " * fieldname * " exists."
                return
            end
            if length(sizedatafield) == 1
                push!(return_outputs[id], clearNP1(fieldname))
                mapping[id][clearNP1(fieldname)] = [fieldname]
            else
                refDof = sizedatafield[2]

                for dof in 1:refDof
                    push!(return_outputs[id], clearNP1(fieldname) * Write_Exodus_Results.get_paraviewCoordinates(dof, refDof))
                    mapping[id][clearNP1(fieldname)*Write_Exodus_Results.get_paraviewCoordinates(dof, refDof)] = fieldname
                end
            end
        end
    end
    return return_outputs, mapping
end

function clearNP1(name)
    if "NP1" == name[end-2:end]
        return name[1:end-3]
    end
    return name
end
end