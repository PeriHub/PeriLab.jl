module IO
include("read_inputdeck.jl")
include("mesh_data.jl")
include("exodus_export.jl")
include("../Support/Parameters/parameter_handling.jl")
using .Read_Input_Deck
using .Read_Mesh
using .Write_Exodus_Results
export close_files
export initialize_data
export init_write_results
export write_results

function close_files(exos)
    @info "Closing output files"
    for exo in exos
        close(exo)
    end
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
                @error "No field " * fieldname * " exists."
                return
            end
            if length(sizedatafield) == 1
                mapping[id][clearNP1(fieldname)] = [fieldname, result_id, 1]
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
    return Read_Mesh.init_data(read_input_file(filename), datamanager, comm)
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
    outputs = get_results_mapping(params, datamanager)
    for i in eachindex(exos)
        exos[i] = Write_Exodus_Results.init_results_in_exodus(exos[i], outputs[i], coords, block_Id, nsets)
    end

    return exos, outputs
end

function read_input_file(filename)
    return Read_Input_Deck.read_input_file(filename)
end

function write_results(exos, step, dt, outputs, datamanager)
    time = (step - 1) * dt

    for id in eachindex(exos)
        exos[id] = write_step_and_time(exos[id], step, time)
        exos[id] = Write_Exodus_Results.write_nodal_results_in_exodus(exo[id], step, outputs[id], datamanager)
    end
end

end