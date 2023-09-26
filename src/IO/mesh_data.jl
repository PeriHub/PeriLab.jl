# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Read_Mesh
using CSV
using DataFrames
using MPI
using TimerOutputs
using LinearAlgebra
using NearestNeighbors: BallTree
using NearestNeighbors: inrange
include("../Support/Parameters/parameter_handling.jl")
include("../MPI_communication/MPI_init.jl")
include("../Support/geometry.jl")
using .Geometry

#export read_mesh
#export load_mesh_and_distribute
export init_data

TOLERANCE = 1.0e-14

function init_data(params, datamanager, comm, to)
    @timeit to "init_data - mesh_data,jl" begin
        ranks = MPI.Comm_size(comm)
        if (MPI.Comm_rank(comm)) == 0
            @timeit to "load_and_evaluate_mesh" distribution, mesh, ntype, overlap_map, nlist, dof = load_and_evaluate_mesh(params, ranks)
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
        @timeit to "get_local_overlap_map" overlap_map = get_local_overlap_map(overlap_map, distribution, ranks)
        @timeit to "distribution_to_cores" datamanager = distribution_to_cores(comm, datamanager, mesh, distribution, dof)
        @timeit to "distribute_neighborhoodlist_to_cores" datamanager = distribute_neighborhoodlist_to_cores(comm, datamanager, nlist, distribution)
        datamanager.set_block_list(datamanager.get_field("Block_Id"))
        datamanager = get_bond_geometry(datamanager) # gives the initial length and bond damage
        @info "Finish init data"
    end
    return datamanager, params
end
"""
    get_local_overlap_map()

    Changes entries in the overlap map from the global numbering to the local computer core one.
    Inputs:
    - `overlap_map::Dict{Int64, Dict{Int64, String}}`: overlap map with global nodes.
    - `distribution Array{Int64}`: global nodes distribution at cores, needed for the gobal to local mapping
    - `ranks Array{Int64}` : number of used computer cores
    Returns:
    - `overlap_map::Dict{Int64, Dict{Int64, String}}`: returns overlap map with local nodes.

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
    return [glob_to_loc[global_node] for global_node in global_nodes if haskey(glob_to_loc, global_node)]
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
    rank = MPI.Comm_rank(comm)
    if rank == 0
        meshdata = check_mesh_elements(mesh, dof)
    else
        meshdata = Dict()
    end
    meshdata = send_value(comm, 0, meshdata)
    for fieldname in keys(meshdata)
        fieldDof = length(meshdata[fieldname]["Mesh ID"])
        datafield = datamanager.create_constant_node_field(fieldname, meshdata[fieldname]["Type"], fieldDof)
        for (localDof, meshID) in enumerate(meshdata[fieldname]["Mesh ID"])
            if rank == 0
                send_msg = meshdata[fieldname]["Type"].(mesh[!, meshID])
                # example send_msg = Float32.(mesh[!, names(mesh)[idof]])
                # redefine from Float64 standard to Float32 for MPI
            else
                send_msg = 0
            end
            if fieldDof == 1
                datafield[:] = send_vector_from_root_to_core_i(comm, send_msg, datafield, distribution)
            else
                datafield[:, localDof] = send_vector_from_root_to_core_i(comm, send_msg, datafield[:, localDof], distribution)
            end
        end
    end
    return datamanager
end
"""
    check_mesh_elements(mesh, dof)

Process and analyze mesh data to create an dictionary containing information
about mesh elements for further processing.

# Arguments
- `mesh::DataFrame`: The input mesh data represented as a DataFrame.
- `dof::Int64`: The degrees of freedom (DOF) for the mesh elements.

# Returns
A dictionary containing information about mesh elements, which can be used for
further processing or uploading.

# Example
```julia
mesh_data = DataFrame(x1 = [1.0, 2.0, 3.0], x2 = [4.0, 5.0, 6.0], volume = [10.0, 20.0, 30.0])
dof = 3
result = check_mesh_elements(mesh_data, dof)
"""
function check_mesh_elements(mesh, dof)
    mnames = names(mesh)
    meshInfoDict = Dict{String,Dict{String,Any}}()

    for (id, mesh_entry) in enumerate(mnames)
        fieldDof = 1
        if ("y" == mesh_entry) | ("z" == mesh_entry) | (mesh_entry[1:end-1] in keys(meshInfoDict))
            continue
        end

        if "x" == mesh_entry
            name = "Coordinates"
            meshID = ["x", "y"]
            if dof == 3
                meshID = ["x", "y", "z"]
            end
        elseif "volume" == mesh_entry
            name = "Volume"
            meshID = ["volume"]
        elseif "block_id" == mesh_entry
            name = "Block_Id"
            meshID = ["block_id"]
        else
            if 'x' == mesh_entry[end] # -> char instead of string
                if id + 1 <= length(mnames)
                    if mnames[id+1][end] == 'y' # -> char instead of string
                        name = mesh_entry[1:end-1]
                        meshID = [name * "x", name * "y"]
                        if dof == 3
                            meshID = [name * "x", name * "y", name * "z"]
                        end
                    end
                end
            else
                name = mesh_entry
                meshID = [name]
            end
        end

        if typeof(mesh[1, meshID[1]]) == Bool
            vartype = Bool
        else
            vartype = typeof(sum(sum(mesh[:, mid])
                                 for mid in meshID))
        end
        if vartype == Float64
            vartype = Float32
        end
        meshInfoDict[name] = Dict{String,Any}("Mesh ID" => meshID, "Type" => vartype)
    end
    if !(haskey(meshInfoDict, "Coordinates"))
        @error "No coordi defined"
    end
    if !(haskey(meshInfoDict, "Block_Id"))
        @error "No blocks defined"
    end
    if !(haskey(meshInfoDict, "Volume"))
        @error "No volumes defined"
    end
    return meshInfoDict
end

function get_header(filename)
    file = open(filename, "r")
    header_line = 0
    for line in eachline(file)#
        header_line += 1
        if contains(line, "header:")
            close(file)
            return header_line, convert(Vector{String}, split(line[9:end], ' '))
        end
    end
end

function read_mesh(filename::String)
    if !isfile(filename)
        @error "File $filename does not exist"
    end
    @info "Read Mesh File $filename"
    header_line, header = get_header(filename)
    return CSV.read(filename, DataFrame; delim=" ", header=header, skipto=header_line + 1)
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
    nlist = apply_bond_filters(nlist, mesh, params, dof)
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


#dict instead of arrays, because it is more readable


function _init_overlap_map_(size)
    #[[[], [[2], [1]], [[2], [5]]], [[[1], [3]], [], [[1, 2], [6, 7]]], [[], [[1, 2], [4, 5]], []]]
    overlap_map = Dict{Int64,Dict{Int64,Dict{String,Vector{Int64}}}}()
    for i in 1:size
        overlap_map[i] = Dict{Int64,Dict{String,Vector{Int64}}}()
        for j in 1:size
            if i != j
                overlap_map[i][j] = Dict{String,Vector{Int64}}("Slave" => Int64[], "Master" => Int64[])
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

function bondIntersectsDisk(p0::Vector{Float64}, p1::Vector{Float64}, center::Vector{Float64}, normal::Vector{Float64}, radius::Float64)
    numerator = sum((center .- p0) .* normal)
    denominator = sum((p1 .- p0) .* normal)
    if abs(denominator) < TOLERANCE
        # Line is parallel to the plane, may or may not lie on the plane
        # If it does lie on the plane, then the numerator will be zero
        # In either case, this function will return "no intersection"
        t = Inf
    else
        # The line intersects the plane
        t = numerator / denominator
    end

    if t < 0.0 || t > 1.0
        return false
    end

    # Intersection point
    x = p0 .+ t .* (p1 .- p0)

    # Check if the intersection point is within the disk
    distance_squared = sum((x .- center) .^ 2)

    if distance_squared < radius^2
        return true
    end

    return false
end

function bondIntersectInfinitePlane(p0::Vector{Float64}, p1::Vector{Float64}, lower_left_corner::Vector{Float64}, normal::Vector{Float64})
    numerator = sum((lower_left_corner .- p0) .* normal)
    denominator = sum((p1 .- p0) .* normal)
    t = 0.0
    if abs(denominator) < TOLERANCE
        # Line is parallel to the plane
        # It may or may not lie on the plane
        # If it does lie on the plane, then the numerator will be zero
        # In either case, this function will return "no intersection"
        t = Inf
    else
        # The line intersects the plane
        t = numerator / denominator
    end

    # Determine if the line segment intersects the plane
    intersection = false
    if 0.0 <= t <= 1.0
        intersection = true
    end

    # Intersection point
    x = p0 .+ t .* (p1 .- p0)

    return intersection, x
end

function bondIntersect(x::Vector{Float64}, lower_left_corner::Vector{Float64}, bottom_unit_vector::Vector{Float64}, normal::Vector{Float64}, side_length::Float64, bottom_length::Float64)
    zero = TOLERANCE
    one = 1.0 + zero
    intersects = false
    ua = cross(bottom_unit_vector, normal)
    r = x
    dr = r .- lower_left_corner
    aa = sum(dr .* ua)
    bb = sum(dr .* bottom_unit_vector)

    if -zero < aa && aa / side_length < one && -zero < bb && bb / bottom_length < one
        intersects = true
    end

    return intersects
end



function apply_bond_filters(nlist, mesh, params, dof)
    bond_filters = get_bond_filters(params)
    if bond_filters[1]
        coor = names(mesh)[1:dof]
        nnodes = length(mesh[!, coor[1]])
        data = zeros(dof, nnodes)
        for i in 1:dof
            data[i, :] = values(mesh[!, coor[i]])
        end
        # balltree = BallTree(data)

        # println(bond_filters)
        for (name, filter) in bond_filters[2]

            if filter["Type"] == "Disk"
                center = [filter["Center_X"], filter["Center_Y"], filter["Center_Z"]]
                normal = [filter["Normal_X"], filter["Normal_Y"], filter["Normal_Z"]]
                for i in 1:nnodes
                    filter_flag = fill(true, length(nlist[i]))
                    for (jId, neighbor) in enumerate(nlist[i])
                        filter_flag[jId] = !bondIntersectsDisk(data[:, i], data[:, neighbor], center, normal, filter["Radius"])
                    end
                    nlist[i] = nlist[i][filter_flag]
                end

            elseif filter["Type"] == "Rectangular_Plane"
                normal = [filter["Normal_X"], filter["Normal_Y"], filter["Normal_Z"]]
                lower_left_corner = [filter["Lower_Left_Corner_X"], filter["Lower_Left_Corner_Y"], filter["Lower_Left_Corner_Z"]]
                bottom_unit_vector = [filter["Bottom_Unit_Vector_X"], filter["Bottom_Unit_Vector_Y"], filter["Bottom_Unit_Vector_Z"]]
                bottom_length = filter["Bottom_Length"]
                side_length = filter["Side_Length"]
                for i in 1:nnodes
                    filter_flag = fill(true, length(nlist[i]))
                    for (jId, neighbor) in enumerate(nlist[i])
                        intersect_inf_plane, x = bondIntersectInfinitePlane(data[:, i], data[:, neighbor], lower_left_corner, normal)
                        bond_intersect = bondIntersect(x, lower_left_corner, bottom_unit_vector, normal, side_length, bottom_length)
                        filter_flag[jId] = !(intersect_inf_plane && bond_intersect)
                    end
                    nlist[i] = nlist[i][filter_flag]
                end
            end
        end
        @info "Finished apply Bond Filters"
    end
    return nlist
end

function glob_to_loc(distribution)
    glob_to_loc = Dict{Int64,Int64}()
    for id in eachindex(distribution)
        glob_to_loc[distribution[id]] = id
    end
    return glob_to_loc
end

end