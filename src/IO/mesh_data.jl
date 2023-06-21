#module mesh_data
using CSV
using DataFrames
using MPI

using NearestNeighbors
include(pwd() * "/Parameter/parameter_handling.jl")
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

        distribute(meshdata, params, ranksize)
    else
        #MPI.Barrier(comm) # notwendig?
    end
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


end

function neighbors(mesh, params, coor)
    nnodes = length(mesh[!, coor[1]])
    dof = length(coor)
    data = zeros(dof, nnodes) # (dof, number of points)
    neighborList = fill([], nnodes)

    for i in 1:dof
        data[i, :] = values(mesh[!, coor[i]]) # memory mapping?
    end
    balltree = BallTree(data)
    for i in 1:nnodes
        neighborList[i] = inrange(balltree, data[:, i], get_horizon(params, values(mesh[!, "block_id"])[i]), true)
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
    MPI.Scatter(Array(data), local_chunk, chunk_size, MPI.INT, MPI.ROOT, MPI.COMM_WORLD)

    # Print the local chunk on each process
    println("Rank $rank: ", local_chunk)
end