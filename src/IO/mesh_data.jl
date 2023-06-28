#module mesh_data
using CSV
using DataFrames
using MPI

using NearestNeighbors
include("../Support/parameter_handling.jl")
#export read_mesh
#export load_mesh_and_distribute
function read_mesh(filename::String)

    header = ["x", "y", "z", "block_id", "volume"]

    if !isfile(filename)
        @error "File $filename does not exist"
    end
    return CSV.read(filename, DataFrame; delim=" ", header=header, skipto=2)
end


function load_mesh_and_distribute(params, comm)

    if (MPI.Comm_rank(comm)) == 0
        meshdata = read_mesh(get_mesh_name(params))
        ranksize = MPI.Comm_size(comm)
        println("$(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")

        topo = distribute(meshdata, params, ranksize)
    else
        #MPI.Barrier(comm) # notwendig?
    end
    return meshdata, topo
end


function distribute(mesh, params, ranksize)
    coor = names(mesh)
    if coor[3] == "z"
        nlist = neighbors(mesh, params, coor[1:3])
    else
        nlist = neighbors(mesh, params, coor[1:2])
    end
    topology = create_topology(nlist, ranksize)
end

function create_topology(nlist, size)

    nnodes = length(nlist)
    if size == 1
        distribution = [collect(1:nnodes)]
    else
        distribution = create_base_chunk(nnodes, size)
        println("distrib $distribution")
        println()
        # check neighborhood & overlap -> all nodes after chunk are overlap
        for i in 1:size
            nchunks = length(distribution[i])
            for j in 1:nchunks
                id = distribution[i, j]
                nneighbors = length(nlist[id])
                for k in 1:nneighbors
                    if length(filter(x -> x == neighbors[nlist[id, k]], distribution[i])) < 1
                        append!(distribution[i], [nlist[id, k]])
                    end
                end
            end
        end

    end
    println("distrib $distribution")
    return distribution
end

function create_base_chunk(nnodes, size)
    # Calculate the initial size of each chunk
    chunk_size = div(nnodes, size)
    # Split the data into chunks
    distribution = fill([], size)
    for i in 1:size
        start_idx = (i - 1) * chunk_size + 1
        end_idx = min(i * chunk_size, nnodes)
        if i == size && end_idx < nnodes
            end_idx = nnodes
        end
        distribution[i] = collect(start_idx:end_idx)
    end
    return distribution
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