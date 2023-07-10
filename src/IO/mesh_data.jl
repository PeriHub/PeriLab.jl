#module mesh_data
using CSV
using DataFrames
using MPI

using NearestNeighbors: BallTree
using NearestNeighbors: inrange
include("../Support/Parameters/parameter_handling.jl")
include("../MPI_communication/MPI_init.jl")

#export read_mesh
#export load_mesh_and_distribute

function init_data(filename, datamanager, comm)

    parameter = read_input_file(filename)

    if (MPI.Comm_rank(comm)) == 0
        distribution, mesh, ntype, overlap_map, dof = load_mesh_and_distribute(parameter, MPI.Comm_size(comm))
        #nodeList = zeros(Int64, nmasters + nslaves)
    else
        nmasters = 0
        nslaves = 0
        overlap_map = nothing
        dof = 0
        mesh = []
        ntype = Dict("masters" => 0, "slaves" => 0)
        distribution = 0
    end
    dof = send_value(comm, 0, dof)
    overlap_map = send_value(comm, 0, overlap_map)
    nmasters::Int64 = send_single_value_from_vector(comm, 0, ntype["masters"], Int64)
    nslaves::Int64 = send_single_value_from_vector(comm, 0, ntype["slaves"], Int64)

    #ntype = send_value(comm, 0, ntype)
    #println(MPI.Comm_rank(comm), " over ", overlap_map, " dof ", dof)
    datamanager.set_nnodes(nmasters + nslaves)
    #println(datamanager.get_nnodes())

    coor = datamanager.create_constant_node_field("Coordinates", Float32, dof)

    for idof in 1:dof
        send_msg = 0

        if MPI.Comm_rank(comm) == 0
            send_msg = mesh[!, names(mesh)[idof]]
        end
        coor[:, idof] = send_vector_from_root_to_core_i(comm, send_msg, coor[:, idof], distribution)
    end

    println(MPI.Comm_rank(comm), " coor ", coor, " dof ", dof)
    data = 0
    return datamanager, parameter
end
function init_distribution_at_cores()
    init_data_field(dof, type)
end
function read_mesh(filename::String)

    header = ["x", "y", "z", "block_id", "volume"]

    if !isfile(filename)
        @error "File $filename does not exist"
    end
    @info "Read File $filename"
    return CSV.read(filename, DataFrame; delim=" ", header=header, skipto=2)
end
function set_dof(mesh)
    if "z" in names(mesh)
        return 3
    else
        return 2
    end
end

function load_mesh_and_distribute(params, ranksize)

    mesh = read_mesh(get_mesh_name(params))
    dof = set_dof(mesh)
    nlist = create_neighborhoodlist(mesh, params, dof)
    distribution, ptc, ntype = node_distribution(nlist, ranksize)

    overlap_map = create_overlap_map(distribution, ptc, ranksize)

    # das kann auf allen Kernen gemacht werden und sollte es auch
    #globToLoc = create_global_to_local_map(distribution)

    # information = Dict("Meshdata" => meshdata, "Nodetype" => ntype, "Overlap_map" => overlap_map, "Node_distribution" => distribution, "Global_to_local" => globToLoc)
    return distribution, mesh, ntype, overlap_map, dof
end

function get_nnodes_per_core(field)
    nnodes = Int64[]
    le = length(field)
    for i in 1:le
        append!(nnodes, field[i])
    end
    return nnodes
end

function create_neighborhoodlist(mesh, params, dof)
    coor = names(mesh)
    nlist = neighbors(mesh, params, coor[1:dof])
    return nlist
end

function node_distribution(nlist, size)

    nnodes = length(nlist)
    ntype = Dict("masters" => Int64[], "slaves" => Int64[])
    if size == 1
        distribution = [collect(1:nnodes)]
        overlap_map = [[[]]]
        ptc = []
        append!(ntype["masters"], nnodes)
        append!(ntype["slaves"], 0)
    else

        distribution, ptc = create_base_chunk(nnodes, size)
        println("distrib $distribution")
        println()
        # check neighborhood & overlap -> all nodes after chunk are overlap
        for i in 1:size
            nchunks = length(distribution[i])
            append!(ntype["masters"], nchunks)

            tempid = Int64[]
            for j in 1:nchunks
                id = distribution[i][j]
                # find all nodes which are not in chunk
                # add them to list as slave nodes for data exchange
                # findall give the index of the valid statement
                # that means that this indices have to be used to obtain the values
                indices = findall(item -> item != i, ptc[nlist[id]])
                if length(indices) > 0
                    append!(tempid, nlist[id][indices])
                end
            end
            # only single new elements where added
            append!(distribution[i], sort(unique(tempid)))
            append!(ntype["slaves"], length(unique(tempid)))
        end

    end
    println("distrib $distribution")
    return distribution, ptc, ntype
end
function create_overlap_map(distribution, ptc, size)
    #[[[], [[2], [1]], [[2], [5]]], [[[1], [3]], [], [[1, 2], [6, 7]]], [[], [[1, 2], [4, 5]], []]]
    overlap_map = _init_overlap_map_(size)
    if size > 1
        for i in 1:size
            vector = distribution[i]
            indices = findall(item -> item != i, ptc[vector])
            le = length(indices)
            from_index = []
            to_index = []
            for j in 1:le
                from_index = findfirst(item -> item == vector[indices[j]], distribution[ptc[vector[indices[j]]]])
                to_index = indices[j]
                append!(overlap_map[ptc[vector[indices[j]]]][i][1], from_index)
                append!(overlap_map[ptc[vector[indices[j]]]][i][2], to_index)
            end

        end
    end
    return overlap_map
end


function _init_overlap_map_(size)
    #[[[], [[2], [1]], [[2], [5]]], [[[1], [3]], [], [[1, 2], [6, 7]]], [[], [[1, 2], [4, 5]], []]]
    overlap_map = []
    for i in 1:size
        append!(overlap_map, [[]])
        for j in 1:size
            append!(overlap_map[i], [[]])
        end
    end
    for i in 1:size
        for j in 1:size
            append!(overlap_map[i][j], [Int64[], Int64[]])
        end
    end
    return overlap_map
end
function create_base_chunk(nnodes, size)
    # Calculate the initial size of each chunk
    # for a nearly equal number of nodes vs. cores this algorithm might lead to the problem, 
    # that the last core is not equally loaded
    chunk_size = div(nnodes, size)
    # Split the data into chunks
    distribution = fill([], size)
    point_to_core = zeros(Int32, nnodes)
    for i in 1:size
        start_idx = (i - 1) * chunk_size + 1
        end_idx = min(i * chunk_size, nnodes)
        if i == size && end_idx < nnodes
            end_idx = nnodes
        end
        distribution[i] = collect(start_idx:end_idx)
        point_to_core[start_idx:end_idx] .= i
    end
    return distribution, point_to_core
end

function neighbors(mesh, params, coor)
    """
    neighbors(mesh, params, coor)

    Compute the neighbor list for each node in a mesh based on their proximity using a BallTree data structure.

    # Arguments
    - `mesh`: A mesh data structure containing the coordinates and other information.
    - `params`: Parameters needed for computing the neighbor list.
    - `coor`: A vector of coordinate names along which to compute the neighbor list.

    # Returns
    An array of neighbor lists, where each element represents the neighbors of a node in the mesh.
    """
    nnodes = length(mesh[!, coor[1]])
    dof = length(coor)
    data = zeros(dof, nnodes)
    neighborList = fill([], nnodes)

    for i in 1:dof
        data[i, :] = values(mesh[!, coor[i]]) # memory mapping?
    end
    balltree = BallTree(data)
    for i in 1:nnodes
        neighborList[i] = inrange(balltree, data[:, i], get_horizon(params, mesh[!, "block_id"][i]), true)
        # avoid self reference in neighborhood
        index = findfirst(x -> x == i, neighborList[i])
        deleteat!(neighborList[i], index)
    end
    return neighborList
end

function init_vectors_for_processes(data, comm, vector)

    # Get the rank and size of the current process
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    size = MPI.Comm_size(MPI.COMM_WORLD)

    # Define the size of the vector
    vector_size = length(data)

    # Calculate the size of each chunk
    chunk_size = vector_size ÷ size

    # Initialize the local chunk of data on each process
    local_chunk = Vector{Int}(undef, chunk_size)

    # Create the complete vector of data on the root process (rank 0)
    if rank == 0
        data = 1:vector_size
    else
        data = Vector{Int}(undef, 0)
    end

    # Scatter the data chunks from the root process to all processes
    #MPI.Scatter(Array(data), local_chunk, chunk_size, MPI.INT, MPI.ROOT, MPI.COMM_WORLD)

    # Print the local chunk on each process
    println("Rank $rank: ", local_chunk)
end