module Read_Mesh
using CSV
using DataFrames
using MPI

using NearestNeighbors: BallTree
using NearestNeighbors: inrange
include("../Support/Parameters/parameter_handling.jl")
include("../MPI_communication/MPI_init.jl")
include("../Support/geometry.jl")
using .Geometry

#export read_mesh
#export load_mesh_and_distribute
export init_data
function init_data(params, datamanager, comm)
    ranks = MPI.Comm_size(comm)
    if (MPI.Comm_rank(comm)) == 0
        distribution, mesh, ntype, overlap_map, nlist, dof = load_and_evaluate_mesh(params, ranks)
    else
        nmasters = 0
        nslaves = 0
        ntype = Dict("masters" => 0, "slaves" => 0)
        nlist = 0
        dof = 0
        mesh = []
        overlap_map = nothing
        distribution = nothing
    end
    dof = send_value(comm, 0, dof)
    dof = datamanager.set_dof(dof)
    overlap_map = send_value(comm, 0, overlap_map)

    distribution = send_value(comm, 0, distribution)
    nlist = send_value(comm, 0, nlist)
    datamanager.set_overlap_map(overlap_map)
    nmasters::Int64 = send_single_value_from_vector(comm, 0, ntype["masters"], Int64)
    nslaves::Int64 = send_single_value_from_vector(comm, 0, ntype["slaves"], Int64)
    datamanager.set_nmasters(nmasters)
    datamanager.set_nslaves(nslaves)
    define_nsets(params, datamanager)
    # defines the order of the global nodes to the local core nodes
    datamanager.set_glob_to_loc(glob_to_loc(distribution[MPI.Comm_rank(comm)+1]))
    overlap_map[:] = get_local_overlap_map(overlap_map, distribution, ranks)
    datamanager = distribution_to_cores(comm, datamanager, mesh, distribution, dof)
    datamanager = distribute_neighborhoodlist_to_cores(comm, datamanager, nlist, distribution)
    datamanager = get_bond_geometry(datamanager) # gives the initial length and bond damage
    @info "Finish init data"
    return datamanager, params
end
"""
    get_local_overlap_map()

    Changes entries in the overlap map from the global numbering to the local computer core one.
    Inputs:
    - `overlap_map` (array): overlap map with global nodes.
    - `distribution` (array): global nodes distribution at cores, needed for the gobal to local mapping
    - `ranks` (Int): number of used computer cores
    Returns:
    - `overlap_map` (array): returns overlap map with local nodes.

    Example:
    ```julia
    get_local_overlap_map(overlap_map, distribution, ranks)  # returns local nodes 
    ```
    """
function get_local_overlap_map(overlap_map, distribution, ranks)
    if ranks == 1
        return overlap_map
    end
    for irank in 1:ranks
        ilocal = glob_to_loc(distribution[irank])
        for jrank in 1:ranks
            if irank != jrank
                overlap_map[irank][jrank]["Slave"][:] = local_nodes_from_dict(ilocal, overlap_map[irank][jrank]["Slave"])
                overlap_map[irank][jrank]["Master"][:] = local_nodes_from_dict(ilocal, overlap_map[irank][jrank]["Master"])
            end
        end
    end
    return sort(overlap_map)
end

function local_nodes_from_dict(glob_to_loc, global_nodes)
    return [glob_to_loc[global_node] for global_node in global_nodes if global_node in keys(glob_to_loc)]
end

function distribute_neighborhoodlist_to_cores(comm, datamanager, nlist, distribution)
    send_msg = 0
    lenNlist = datamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nnodes = datamanager.get_nnodes() # master nodes
    rank = MPI.Comm_rank(comm)
    if rank == 0
        send_msg = get_number_of_neighbornodes(nlist)
    end
    lenNlist[:] = send_vector_from_root_to_core_i(comm, send_msg, lenNlist, distribution)
    nlistCore = datamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    # provide neighborhood only for master nodes
    nlistCore[:] = nlist[distribution[rank+1][1:nnodes]]
    nlistCore[:] = get_local_neighbors(datamanager.get_local_nodes, nlistCore)
    nlist = 0
    return datamanager
end

function get_local_neighbors(mapping, nlistCore)
    for id in eachindex(nlistCore)
        nlistCore[id] = mapping(nlistCore[id])
    end
    return nlistCore
end

function get_bond_geometry(datamanager)
    dof = datamanager.get_dof()
    nnodes = datamanager.get_nnodes()
    nlist = datamanager.get_field("Neighborhoodlist")
    coor = datamanager.get_field("Coordinates")
    bondgeom = datamanager.create_constant_bond_field("Bond Geometry", Float32, dof + 1)
    bondDamage = datamanager.create_constant_bond_field("Bond Damage", Float32, 1)
    bondgeom = Geometry.bond_geometry(1:nnodes, dof, nlist, coor, bondgeom)
    return datamanager
end

function define_nsets(params, datamanager)
    nsets = get_node_sets(params)
    for nset in keys(nsets)
        datamanager.set_nset(nset, nsets[nset])
    end
end

function distribution_to_cores(comm, datamanager, mesh, distribution, dof)
    # init blockID field

    blockID = datamanager.create_constant_node_field("Block_Id", Int64, 1)
    # set value for all cores as send_msg
    # only used for rank = 0
    send_msg = 0
    # distribute blocks
    if MPI.Comm_rank(comm) == 0
        mnames = names(mesh)
        if ("block_id" in mnames) == false
            @error "No blocks defined"
        end
        send_msg = mesh[!, "block_id"]
    end
    # must be [:] -> to map it in datamanager
    blockID[:] = send_vector_from_root_to_core_i(comm, send_msg, blockID, distribution)
    datamanager.set_block_list(blockID)
    # init coordinate field
    coor = datamanager.create_constant_node_field("Coordinates", Float32, dof)
    # distribute coordinates
    for idof in 1:dof
        if MPI.Comm_rank(comm) == 0
            send_msg = Float32.(mesh[!, names(mesh)[idof]])
        end
        coor[:, idof] = send_vector_from_root_to_core_i(comm, send_msg, coor[:, idof], distribution)
    end
    volume = datamanager.create_constant_node_field("Volume", Float32, 1)
    # distribute blocks
    if MPI.Comm_rank(comm) == 0
        send_msg = Float32.(mesh[!, "volume"]) # because mesh is read as Float64
    end
    volume[:] = send_vector_from_root_to_core_i(comm, send_msg, volume, distribution)
    # additional fields as angles and pointtime can be add automatically

    # send header to all cores
    # create fields
    # check in advance x,yz for dof>1 fields
    # transfer data
    # upload=Dict()
    # if MPI.Comm_rank(comm) == 0
    #     upload=check_elements(mesh, dof)
    # end
    # upload = send_vectors_to_cores(comm, 0, upload)
    # #for name in mesh[!, "volume"]
    # #    send_msg = Float32.(mesh[!, "volume"])
    # end
    # for upload in uploadDict
    #     var = datamanager.create_constant_node_field(upload["Name"], upload["Type"], upload["Dof"])
    # end
    return datamanager
end

function check_elements(mesh, dof)
    mnames = names(mesh)

    uploadDict = []
    count = 0
    for nam in mnames
        count += 1
        fieldDof = 1
        if "x" == nam
            name = "Coordinates"
            meshentry = collect(count:count+dof)
            fieldDof = dof

        elseif "volume" == nam
            name = "Volume"
            meshentry = collect(count:count)

        elseif "block_id" == nam
            name = "Block_Id"
            meshentry = collect(count:count)
        else
            if "x" in nam
                meshentry = collect(count:count+dof)
                name = nam[1:length(nam)-2]
            else
                meshentry = collect(count:count)
                name = nam
            end
        end

        append!(uploadDict, Dict("Name" => name, "Pos" => meshentry, "Type" => typeof(mesh[nam][1]), "Dof" => fieldDof))

    end
    return uploadDict
end

function get_header(filename)
    file = open(filename, "r")
    for line in eachline(file)
        if contains(line, '#')
            close(file)
            return convert(Vector{String}, split(line[3:end], ' '))
        end
    end
end

function read_mesh(filename::String)
    if !isfile(filename)
        @error "File $filename does not exist"
    end
    @info "Read Mesh File $filename"
    header = get_header(filename)
    return CSV.read(filename, DataFrame; delim=" ", header=header, skipto=2)
end

function set_dof(mesh)
    if "z" in names(mesh)
        return 3
    else
        @info "2d problem assumed, please define plane stress or plane strain if needed in Physics"
        return 2
    end
end

function load_and_evaluate_mesh(params, ranksize)

    mesh = read_mesh(get_mesh_name(params))
    dof = set_dof(mesh)
    nlist = create_neighborhoodlist(mesh, params, dof)
    @info "Start distribution"
    distribution, ptc, ntype = node_distribution(nlist, ranksize)
    @info "Finished distribution"
    @info "Create Overlap"
    overlap_map = create_overlap_map(distribution, ptc, ranksize)
    @info "Finished Overlap"
    # das kann auf allen Kernen gemacht werden und sollte es auch
    #globToLoc = create_global_to_local_map(distribution)

    # information = Dict("Meshdata" => meshdata, "Nodetype" => ntype, "Overlap_map" => overlap_map, "Node_distribution" => distribution, "Global_to_local" => globToLoc)
    return distribution, mesh, ntype, overlap_map, nlist, dof
end

function get_nnodes_per_core(field)
    nnodes = Int64[]
    for i in 1:eachindex(field)
        append!(nnodes, field[i])
    end
    return nnodes
end

function create_neighborhoodlist(mesh, params, dof)
    coor = names(mesh)
    nlist = neighbors(mesh, params, coor[1:dof])
    @info "Finished init Neighborhoodlist"
    return nlist
end

function get_number_of_neighbornodes(nlist)
    len = length(nlist)
    lenNlist = zeros(Int64, len)
    for id in 1:len
        lenNlist[id] = length(nlist[id])
    end
    return lenNlist
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
    return distribution, ptc, ntype
end

"""
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
"""


#dict instead of arrays, because it is more readable


function _init_overlap_map_(size)
    #[[[], [[2], [1]], [[2], [5]]], [[[1], [3]], [], [[1, 2], [6, 7]]], [[], [[1, 2], [4, 5]], []]]
    overlap_map = Dict{Int64,Dict{Int64,Dict{String,Vector{Int64}}}}()
    for i in 1:size
        overlap_map[i] = Dict{Int64,Dict{String,Vector{Int64}}}()
        for j in 1:size
            if i != j
                overlap_map[i][j] = Dict{String,Vector{Int64}}("Send" => Int64[], "Receive" => Int64[])
            end
        end
    end

    return overlap_map
end
"""
    ptc - point to core map; it gives the cores where the master nodes are -> basis chunk gives this directly


"""
function create_overlap_map(distribution, ptc, size)

    overlap_map = _init_overlap_map_(size)
    if size == 1
        return overlap_map
    end
    for icoreID in 1:size
        # distribution of nodes at core i
        vector = distribution[icoreID]
        # gives core ids of all nodes not master at core icoreID
        for jcoreID in 1:size
            if icoreID == jcoreID
                continue
            end
            indices = findall(item -> item == jcoreID, ptc[vector])
            overlap_map[icoreID][jcoreID]["Slave"] = vector[indices]
            overlap_map[jcoreID][icoreID]["Master"] = vector[indices]
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
    - `params`: paramss needed for computing the neighbor list.
    - `coor`: A vector of coordinate names along which to compute the neighbor list.

    # Returns
    An array of neighbor lists, where each element represents the neighbors of a node in the mesh.
    """
    @info "Init Neighborhoodlist"
    nnodes = length(mesh[!, coor[1]])
    dof = length(coor)
    data = zeros(dof, nnodes)
    neighborList = fill([], nnodes)

    for i in 1:dof
        data[i, :] = values(mesh[!, coor[i]])
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

function glob_to_loc(distribution)
    glob_to_loc = Dict{Int64,Int64}()
    for id in eachindex(distribution)
        glob_to_loc[distribution[id]] = id
    end
    return glob_to_loc
end

end