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
include("../Support/helpers.jl")
using .Geometry

#export read_mesh
#export load_mesh_and_distribute
export init_data

const TOLERANCE = 1.0e-14

"""
    init_data(params::Dict, path::String, datamanager::Module, comm::MPI.Comm, to::TimerOutput)

Initializes the data for the mesh.

# Arguments
- `params::Dict`: The parameters for the simulation.
- `path::String`: The path to the mesh file.
- `datamanager::Data_manager`: The data manager.
- `comm::MPI.Comm`: The MPI communicator.
- `to::TimerOutput`: The timer output.
# Returns
- `datamanager::Data_manager`: The data manager.
- `params::Dict`: The parameters for the simulation.
"""
function init_data(params::Dict, path::String, datamanager::Module, comm::MPI.Comm, to::TimerOutput)
    @timeit to "init_data - mesh_data,jl" begin
        ranks = MPI.Comm_size(comm)
        if (MPI.Comm_rank(comm)) == 0
            @timeit to "load_and_evaluate_mesh" distribution, mesh, ntype, overlap_map, nlist, dof = load_and_evaluate_mesh(params::Dict, path, ranks)
        else
            num_controller = 0
            num_responder = 0
            ntype = Dict("controllers" => 0, "responder" => 0)
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
        num_controller::Int64 = send_single_value_from_vector(comm, 0, ntype["controllers"], Int64)
        num_responder::Int64 = send_single_value_from_vector(comm, 0, ntype["responder"], Int64)
        datamanager.set_num_controller(num_controller)
        datamanager.set_num_responder(num_responder)
        @info "Get node sets"
        define_nsets(params, path, datamanager)
        # defines the order of the global nodes to the local core nodes
        datamanager.set_distribution(distribution[MPI.Comm_rank(comm)+1])
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

# Arguments
- `overlap_map::Dict{Int64, Dict{Int64, String}}`: overlap map with global nodes.
- `distribution::Vector{Vector{Int64}}`: global nodes distribution at cores, needed for the gobal to local mapping
- `ranks Array{Int64}` : number of used computer cores
# Returns
- `overlap_map::Dict{Int64, Dict{Int64, String}}`: returns overlap map with local nodes.

Example:
```julia
get_local_overlap_map(overlap_map, distribution, ranks)  # returns local nodes 
```
"""
function get_local_overlap_map(overlap_map, distribution::Vector{Vector{Int64}}, ranks::Int64)
    if ranks == 1
        return overlap_map
    end
    for irank in 1:ranks
        ilocal = glob_to_loc(distribution[irank])
        for jrank in 1:ranks
            if irank != jrank
                overlap_map[irank][jrank]["Responder"][:] = local_nodes_from_dict(ilocal, overlap_map[irank][jrank]["Responder"])
                overlap_map[irank][jrank]["Controller"][:] = local_nodes_from_dict(ilocal, overlap_map[irank][jrank]["Controller"])
            end
        end
    end
    return sort(overlap_map)
end

"""
    local_nodes_from_dict(glob_to_loc::Dict{Int,Int}, global_nodes::Vector{Int64})

Changes entries in the overlap map from the global numbering to the local computer core one.

# Arguments
- `glob_to_loc::Dict{Int,Int}`: global to local mapping
- `global_nodes::Vector{Int64}`: global nodes
# Returns
- `overlap_map::Dict{Int64, Dict{Int64, String}}`: returns overlap map with local nodes.
"""
function local_nodes_from_dict(glob_to_loc::Dict{Int,Int}, global_nodes::Vector{Int64})
    return Int64[glob_to_loc[global_node] for global_node in global_nodes if haskey(glob_to_loc, global_node)]
end

"""
    distribute_neighborhoodlist_to_cores(comm::MPI.Comm, datamanager::Module, nlist, distribution)

Distributes the neighborhood list to the cores.

# Arguments
- `comm::MPI.Comm`: MPI communicator
- `datamanager::Module`: Data manager
- `nlist`: neighborhood list
- `distribution Array{Int64}`: global nodes distribution at cores
# Returns
- `datamanager::Module`: data manager
"""
function distribute_neighborhoodlist_to_cores(comm::MPI.Comm, datamanager::Module, nlist::Vector{Vector{Int64}}, distribution::Vector{Vector{Int64}})
    send_msg = 0
    lenNlist = datamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    rank = MPI.Comm_rank(comm)
    if rank == 0
        send_msg = get_number_of_neighbornodes(nlist)
    end
    lenNlist[:] = send_vector_from_root_to_core_i(comm, send_msg, lenNlist, distribution)
    nlistCore = datamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)

    nlistCore[:] = nlist[distribution[rank+1][:]]
    nlistCore[:] = get_local_neighbors(datamanager.get_local_nodes, nlistCore)
    nlist = 0
    return datamanager
end

"""
    get_local_neighbors(mapping, nlistCore)

Gets the local neighborhood list from the global neighborhood list

# Arguments
- `mapping`: mapping function
- `nlistCore`: global neighborhood list
# Returns
- `nlistCore`: local neighborhood list
"""
function get_local_neighbors(mapping, nlistCore)
    for id in eachindex(nlistCore)
        nlistCore[id] = mapping(nlistCore[id])
    end
    return nlistCore
end

"""
    get_bond_geometry(datamanager::Module)

Gets the bond geometry

# Arguments
- `datamanager::Module`: Data manager
# Returns
- `datamanager::Module`: data manager
"""
function get_bond_geometry(datamanager::Module)
    dof = datamanager.get_dof()
    nnodes = datamanager.get_nnodes()
    nlist = datamanager.get_field("Neighborhoodlist")
    coor = datamanager.get_field("Coordinates")
    undeformed_bond = datamanager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    bond_damage = datamanager.create_constant_bond_field("Bond Damage", Float64, 1)
    undeformed_bond = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, undeformed_bond)
    return datamanager
end

"""
    define_nsets(params::Dict, path::String, datamanager::Module)

Defines the node sets

# Arguments
- `params::Dict`: Parameters
- `path::String`: Path
- `datamanager::Module`: Data manager
"""
function define_nsets(params::Dict, path::String, datamanager::Module)
    nsets = get_node_sets(params, path)
    for nset in keys(nsets)
        datamanager.set_nset(nset, nsets[nset])
    end
end

"""
    distribution_to_cores(comm::MPI.Comm, datamanager::Module, mesh, distribution, dof::Int64)

Distributes the mesh data to the cores

# Arguments
- `comm::MPI.Comm`: MPI communicator
- `datamanager::Module`: Data manager
- `mesh`: Mesh
- `distribution Array{Int64}`: global nodes distribution at cores
- `dof::Int64`: Degree of freedom
# Returns
- `datamanager::Module`: data manager
"""
function distribution_to_cores(comm::MPI.Comm, datamanager::Module, mesh::DataFrame, distribution::Vector{Vector{Int64}}, dof::Int64)
    # init block_id field
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
                # example send_msg = Float64.(mesh[!, names(mesh)[idof]])
                # redefine from Float64 standard to Float64 for MPI
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
function check_mesh_elements(mesh::DataFrame, dof::Int64)
    mnames = names(mesh)
    meshInfoDict = Dict{String,Dict{String,Any}}()

    for (id, mesh_entry) in enumerate(mnames)
        fieldDof = 1
        if ("y" == mesh_entry) || ("z" == mesh_entry) || (mesh_entry[1:end-1] in keys(meshInfoDict))
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
            if "_x" == mesh_entry[end-1:end]
                if id + 1 <= length(mnames)
                    if mnames[id+1][end-1:end] == "_y"
                        name = mesh_entry[1:end-2]
                        meshID = [name * "_x", name * "_y"]
                        if dof == 3
                            meshID = [name * "_x", name * "_y", name * "_z"]
                        end
                    end
                end
            else
                name = mesh_entry
                meshID = [name]
            end
        end

        if mesh[1, meshID[1]] isa Bool
            vartype = Bool
        else
            vartype = typeof(sum(sum(mesh[:, mid])
                                 for mid in meshID))
        end
        meshInfoDict[name] = Dict{String,Any}("Mesh ID" => meshID, "Type" => vartype)
    end
    if !(haskey(meshInfoDict, "Coordinates"))
        @error "No coordinates defined"
    end
    if !(haskey(meshInfoDict, "Block_Id"))
        @error "No blocks defined"
    end
    if !(haskey(meshInfoDict, "Volume"))
        @error "No volumes defined"
    end
    return meshInfoDict
end

"""
    read_external_topology(filename::String)

Read external topoloy data from a file and return it as a DataFrame.

# Arguments
- `filename::String`: The path to the mesh file.
# Returns
- `external_topology::DataFrame`: The external topology data as a DataFrame.
"""
function read_external_topology(filename::String)
    if !isfile(filename)
        return nothing
    end
    @info "Read external topology file $filename"
    header_line, header = get_header(filename)
    return CSV.read(filename, DataFrame; delim=" ", ignorerepeated=true, header=header, skipto=header_line + 1, comment="#")
end

"""
    read_mesh(filename::String)

Read mesh data from a file and return it as a DataFrame.

# Arguments
- `filename::String`: The path to the mesh file.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
"""
function read_mesh(filename::String)
    if !isfile(filename)
        @error "File $filename does not exist"
        return nothing
    end
    @info "Read mesh file $filename"
    header_line, header = get_header(filename)
    return CSV.read(filename, DataFrame; delim=" ", ignorerepeated=true, header=header, skipto=header_line + 1, comment="#")
end

"""
    set_dof(mesh::DataFrame)

Set the degrees of freedom (DOF) for the mesh elements.

# Arguments
- `mesh::DataFrame`: The input mesh data represented as a DataFrame.
# Returns
- `dof::Int64`: The degrees of freedom (DOF) for the mesh elements.
"""
function set_dof(mesh::DataFrame)
    if "z" in names(mesh)
        return Int64(3)
    end
    @info "2d problem assumed, please define plane stress or plane strain if needed in Physics"
    return Int64(2)
end

"""
    load_and_evaluate_mesh(params::Dict, path::String, ranksize::Int64)

Load and evaluate the mesh data.

# Arguments
- `params::Dict`: The input parameters.
- `path::String`: The path to the mesh file.
- `ranksize::Int64`: The number of ranks.
# Returns
- `distribution::Array{Int64,1}`: The distribution of the mesh elements.
- `mesh::DataFrame`: The mesh data as a DataFrame.
- `ntype::Dict`: The type of the mesh elements.
- `overlap_map::Array{Array{Int64,1},1}`: The overlap map of the mesh elements.
- `nlist::Array{Array{Int64,1},1}`: The neighborhood list of the mesh elements.
- `dof::Int64`: The degrees of freedom (DOF) for the mesh elements.
"""
function load_and_evaluate_mesh(params::Dict, path::String, ranksize::Int64)

    mesh = read_mesh(joinpath(path, get_mesh_name(params)))
    duplicates = findall(nonunique(mesh))
    if length(duplicates) > 0
        @error "Mesh contains duplicate nodes! Nodes: $duplicates"
        return nothing
    end
    external_topology = nothing
    if !isnothing(get_external_topology_name(params))
        external_topology = read_external_topology(joinpath(path, get_external_topology_name(params)))
    end
    if !isnothing(external_topology)
        @info "External topology files was read."
    end
    dof::Int64 = set_dof(mesh)
    nlist = create_neighborhoodlist(mesh, params, dof)
    nlist = apply_bond_filters(nlist, mesh, params, dof)
    if !isnothing(external_topology)
        @info "Create a consistent neighborhood list with external topology definition."
        nlist, topology = create_consistent_neighborhoodlist(external_topology, params["Discretization"]["Input External Topology"], nlist, dof)
    end
    @info "Start distribution"
    distribution, ptc, ntype = node_distribution(nlist, ranksize)
    if haskey(params, "FEM") && !isnothing(external_topology)
        element_distribution, pte = element_distribution(topology, ptc, ranksize)
    end
    @info "Finished distribution"
    @info "Create Overlap"
    overlap_map = create_overlap_map(distribution, ptc, ranksize)
    @info "Finished Overlap"
    @info "Mesh input overview"
    @info "-------------------"
    @info "Number of nodes: $(length(mesh[!, "x"]))"
    @info "Geometrical degrees of freedoms: $dof"
    @info "-------------------"
    return distribution, mesh, ntype, overlap_map, nlist, dof
end

function create_consistent_neighborhoodlist(external_topology::DataFrame, params::Dict, nlist::Vector{Vector{Int64}}, dof::Int64)
    pd_neighbors::Bool = false
    if haskey(params, "Add Neighbor Search")
        pd_neighbors = params["Add Neighbor Search"]
    end
    number_of_elements = length(external_topology[:, 1])
    topology::Vector{Vector{Int64}} = []
    for i_el in 1:number_of_elements
        push!(topology, collect(skipmissing(external_topology[i_el, :])))
    end
    nodes_to_element = [Any[] for _ in 1:maximum(maximum(topology))]
    fe_nodes = Vector{Int64}()
    for (el_id, topo) in enumerate(topology)
        for node in topo
            push!(nodes_to_element[node], el_id)
            push!(fe_nodes, node)
        end
    end
    fe_nodes = unique(fe_nodes)
    for fe_node in fe_nodes
        if !pd_neighbors
            nlist[fe_node] = Int64[]
        end
        for el_id in nodes_to_element[fe_node]
            append!(nlist[fe_node], topology[el_id])
        end
        nlist[fe_node] = unique(nlist[fe_node])
        nlist[fe_node] = filter(x -> x != fe_node, nlist[fe_node])
    end
    return nlist, topology, nodes_to_element
end

"""
    create_neighborhoodlist(mesh::DataFrame, params::Dict, dof::Int64)

Create the neighborhood list of the mesh elements.

# Arguments
- `mesh::DataFrame`: The input mesh data represented as a DataFrame.
- `params::Dict`: The input parameters.
- `dof::Int64`: The degrees of freedom (DOF) for the mesh elements.
# Returns
- `nlist::Array{Array{Int64,1},1}`: The neighborhood list of the mesh elements.
"""
function create_neighborhoodlist(mesh::DataFrame, params::Dict, dof::Int64)
    coor = names(mesh)
    nlist::Vector{Vector{Int64}} = neighbors(mesh, params, coor[1:dof])
    @info "Finished init Neighborhoodlist"
    return nlist
end

"""
    get_number_of_neighbornodes(nlist::Vector{Vector{Int64}})

Get the number of neighbors for each node.

# Arguments
- `nlist::Vector{Vector{Int64}}`: The neighborhood list of the mesh elements.
# Returns
- `lenNlist::Vector{Int64}`: The number of neighbors for each node.
"""
function get_number_of_neighbornodes(nlist::Vector{Vector{Int64}})
    len = length(nlist)
    lenNlist = zeros(Int64, len)
    for id in 1:len
        if length(nlist[id]) == 0
            @error "Node $id has no neighbors please check the horizon."
            return nothing
        end
        lenNlist[id] = length(nlist[id])
    end
    return lenNlist
end

"""
    element_distribution(topology::Vector{Vector{Int64}}, ptc::Vector{Int64}, size::Int64)

Create the distribution of the finite elements. Is needed to avoid multiple element calls. Each element should run only one time at the cores.

# Arguments
- `topology::Vector{Vector{Int64}}`: The topology list of the mesh elements.
- `nlist::Vector{Vector{Int64}}`: The neighborhood list of the mesh elements.
- `size::Int64`: The number of ranks.
# Returns
- `distribution::Vector{Vector{Int64}}`: The distribution of the nodes.
- `etc::Vector{Int64}`: The number of nodes in each rank.
"""
function element_distribution(topology::Vector{Vector{Int64}}, ptc::Vector{Int64}, size::Int64)
    nelements = length(topology)
    if size == 1
        distribution = [collect(1:nelements)]
        etc::Vector{Int64} = []
    else
        distribution, etc = create_base_chunk(nelements, size)
        # check if at least one node of an element is at the same core
        temp = []
        for i_core in 1:size
            push!(temp, Vector{Int64}([]))
            #nchunks = length(distribution[i])
            for el_id in distribution[i_core]
                if !(i_core in ptc[topology[el_id]])
                    push!(temp[i_core], el_id)
                end
            end
        end
        for i_core in 1:size
            if length(temp[i_core]) > 0
                for el_id in temp[i_core]
                    # find core with the lowest number of elements on it
                    # first find cores where all the nodes are
                    # check the number of elements at all of these cores
                    # find the core with the lowest number of elements
                    # put the element there
                    min_core = ptc[topology[el_id][argmin(length.(distribution[ptc[topology[el_id]]]))]]
                    push!(distribution[min_core], el_id)
                    distribution[etc[el_id]] = filter(x -> x != el_id, distribution[etc[el_id]])
                    etc[el_id] = min_core
                end
            end
        end
    end
    return distribution, etc
end

"""
    node_distribution(nlist::Vector{Vector{Int64}}, size::Int64)

Create the distribution of the nodes.

# Arguments
- `nlist::Vector{Vector{Int64}}`: The neighborhood list of the mesh elements.
- `size::Int64`: The number of ranks.
# Returns
- `distribution::Vector{Vector{Int64}}`: The distribution of the nodes.
- `ptc::Vector{Int64}`: Defines at which core / rank each node lies.
- `ntype::Dict`: The type of the nodes.
"""
function node_distribution(nlist::Vector{Vector{Int64}}, size::Int64)

    nnodes = length(nlist)
    ntype = Dict("controllers" => Int64[], "responder" => Int64[])
    if size == 1
        distribution = [collect(1:nnodes)]
        overlap_map = [[[]]]
        ptc::Vector{Int64} = []
        append!(ntype["controllers"], nnodes)
        append!(ntype["responder"], 0)
    else

        distribution, ptc = create_base_chunk(nnodes, size)

        # check neighborhood & overlap -> all nodes after chunk are overlap
        for i in 1:size
            nchunks = length(distribution[i])
            append!(ntype["controllers"], nchunks)

            tempid = Int64[]
            for j in 1:nchunks
                id = distribution[i][j]
                # find all nodes which are not in chunk
                # add them to list as responder nodes for data exchange
                # findall give the index of the valid statement
                # that means that this indices have to be used to obtain the values
                indices = findall(item -> item != i, ptc[nlist[id]])
                if length(indices) > 0
                    append!(tempid, nlist[id][indices])
                end
            end
            # only single new elements where added
            append!(distribution[i], sort(unique(tempid)))
            append!(ntype["responder"], length(unique(tempid)))
        end

    end
    return distribution, ptc, ntype
end


#dict instead of arrays, because it is more readable


"""
    _init_overlap_map_(size)

Initialize the overlap map.

# Arguments
- `size::Int64`: The number of ranks.
# Returns
- `overlap_map::Dict{Int64,Dict{Int64,Dict{String,Vector{Int64}}}}`: The overlap map.
"""
function _init_overlap_map_(size)
    #[[[], [[2], [1]], [[2], [5]]], [[[1], [3]], [], [[1, 2], [6, 7]]], [[], [[1, 2], [4, 5]], []]]
    overlap_map = Dict{Int64,Dict{Int64,Dict{String,Vector{Int64}}}}()
    for i in 1:size
        overlap_map[i] = Dict{Int64,Dict{String,Vector{Int64}}}()
        for j in 1:size
            if i != j
                overlap_map[i][j] = Dict{String,Vector{Int64}}("Responder" => Int64[], "Controller" => Int64[])
            end
        end
    end

    return overlap_map
end

"""
    create_overlap_map(distribution, ptc, size)

Create the overlap map.

# Arguments
- `distribution::Array{Int64,1}`: The distribution of the nodes.
- `ptc::Array{Int64,1}`: The number of nodes in each rank.
- `size::Int64`: The number of ranks.
# Returns
- `overlap_map::Dict{Int64,Dict{Int64,Dict{String,Vector{Int64}}}}`: The overlap map.
"""
function create_overlap_map(distribution::Vector{Vector{Int64}}, ptc::Vector{Int64}, size::Int64)

    overlap_map = _init_overlap_map_(size)
    if size == 1
        return overlap_map
    end
    for icoreID in 1:size
        # distribution of nodes at core i
        vector = distribution[icoreID]
        # gives core ids of all nodes not controller at core icoreID
        for jcoreID in 1:size
            if icoreID == jcoreID
                continue
            end
            indices = findall(item -> item == jcoreID, ptc[vector])
            overlap_map[icoreID][jcoreID]["Responder"] = vector[indices]
            overlap_map[jcoreID][icoreID]["Controller"] = vector[indices]
        end
    end
    return overlap_map
end

"""
    create_base_chunk(nnodes::Int64, size::Int64)

Calculate the initial size of each chunk for a nearly equal number of nodes vs. cores this algorithm might lead to the problem, that the last core is not equally loaded

# Arguments
- `nnodes::Int64`: The number of nodes.
- `size::Int64`: The number of cores.
# Returns
- `distribution::Array{Int64,1}`: The distribution of the nodes.
- `point_to_core::Array{Int64,1}`: The number of nodes in each rank.
"""
function create_base_chunk(nnodes::Int64, size::Int64)
    if size > nnodes
        @error "Number of cores $size exceeds number of nodes $nnodes."
        return nothing, nothing
    end
    chunk_size = div(nnodes, size)
    # Split the data into chunks
    distribution = fill(Int64[], size)
    point_to_core::Vector{Int64} = zeros(Int64, nnodes)
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

"""
    neighbors(mesh, params::Dict, coor)

Compute the neighbor list for each node in a mesh based on their proximity using a BallTree data structure.

# Arguments
- `mesh`: A mesh data structure containing the coordinates and other information.
- `params`: paramss needed for computing the neighbor list.
- `coor`: A vector of coordinate names along which to compute the neighbor list.

# Returns
An array of neighbor lists, where each element represents the neighbors of a node in the mesh.
"""
function neighbors(mesh::DataFrame, params::Dict, coor::Union{Vector{Int64},Vector{String}})
    @info "Init Neighborhoodlist"
    nnodes = length(mesh[!, coor[1]])
    dof = length(coor)
    data = zeros(dof, nnodes)
    neighborList = fill(Vector{Int64}([]), nnodes)

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

"""
    bondIntersectsDisk(p0::Vector{Float64}, p1::Vector{Float64}, center::Vector{Float64}, normal::Vector{Float64}, radius::Float64)

Check if a line segment intersects a disk.

# Arguments
- `p0::Vector{Float64}`: The start point of the line segment.
- `p1::Vector{Float64}`: The end point of the line segment.
- `center::Vector{Float64}`: The center of the disk.
- `normal::Vector{Float64}`: The normal of the plane.
- `radius::Float64`: The radius of the disk.
# Returns
- `Bool`: True if the line segment intersects the disk, False otherwise.
"""
function bondIntersectsDisk(p0::Vector{Float64}, p1::Vector{Float64}, center::Vector{Float64}, normal::Vector{Float64}, radius::Float64)
    numerator = dot((lower_left_corner - p0), normal)
    denominator = dot((p1 - p0), normal)
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
    x = p0 + t .* (p1 - p0)

    # Check if the intersection point is within the disk
    distance = norm(x - center)

    if abs(distance) < radius
        return true
    end

    return false
end

"""
    bondIntersectInfinitePlane(p0::Vector{Float64}, p1::Vector{Float64}, lower_left_corner::Vector{Float64}, normal::Vector{Float64})

Check if a line segment intersects an infinite plane.

# Arguments
- `p0::Vector{Float64}`: The start point of the line segment.
- `p1::Vector{Float64}`: The end point of the line segment.
- `lower_left_corner::Vector{Float64}`: The lower left corner of the plane.
- `normal::Vector{Float64}`: The normal of the plane.
# Returns
- `Bool`: True if the line segment intersects the plane, False otherwise.
"""
function bondIntersectInfinitePlane(p0::Vector{Float64}, p1::Vector{Float64}, lower_left_corner::Vector{Float64}, normal::Vector{Float64})
    numerator = dot((lower_left_corner - p0), normal)
    denominator = dot((p1 - p0), normal)
    if abs(denominator) < TOLERANCE
        # Line is parallel to the plane
        # It may or may not lie on the plane
        # If it does lie on the plane, then the numerator will be zero
        # In either case, this function will return "no intersection"
        return false, undef
    end
    # The line intersects the plane
    t = numerator / denominator

    # Determine if the line segment intersects the plane
    if 0.0 <= t <= 1.0
        return true, p0 + t .* (p1 - p0)
    end
    # Intersection point
    return false, undef
end

"""
    bondIntersectRectanglePlane(x::Vector{Float64}, lower_left_corner::Vector{Float64}, bottom_unit_vector::Vector{Float64}, normal::Vector{Float64}, side_length::Float64, bottom_length::Float64)

Check if a bond intersects a rectangle plane.

# Arguments
- `x::Vector{Float64}`: The point.
- `lower_left_corner::Vector{Float64}`: The lower left corner of the rectangle.
- `bottom_unit_vector::Vector{Float64}`: The unit vector along the bottom of the rectangle.
- `normal::Vector{Float64}`: The normal of the plane.
- `side_length::Float64`: The side length of the rectangle.
- `bottom_length::Float64`: The bottom length of the rectangle.
# Returns
- `Bool`: True if the point is inside the rectangle, False otherwise.
"""
function bondIntersectRectanglePlane(x::Vector{Float64}, lower_left_corner::Vector{Float64}, bottom_unit_vector::Vector{Float64}, normal::Vector{Float64}, side_length::Float64, bottom_length::Float64)
    zero = TOLERANCE
    one = 1.0 + zero

    dr = x - lower_left_corner
    bb = dot(dr, bottom_unit_vector)
    if -zero < bb && bb / bottom_length < one
        if length(normal) == 2
            return true
        end
        ua = cross(bottom_unit_vector, normal)
        aa = dot(dr, ua)

        if -zero < aa && aa / side_length < one
            return true
        end
    end

    return false
end

"""
    apply_bond_filters(nlist::Vector{Vector{Int64}}, mesh::DataFrame, params::Dict, dof::Int64)

Apply the bond filters to the neighborhood list.

# Arguments
- `nlist::Vector{Vector{Int64}}`: The neighborhood list.
- `mesh::DataFrame`: The mesh.
- `params::Dict`: The parameters.
- `dof::Int64`: The degrees of freedom.
# Returns
- `nlist::Vector{Vector{Int64}}`: The filtered neighborhood list.
"""
function apply_bond_filters(nlist::Vector{Vector{Int64}}, mesh::DataFrame, params::Dict, dof::Int64)
    bond_filters = get_bond_filters(params)
    if bond_filters[1]
        @info "Apply bond filters"
        coor = names(mesh)[1:dof]
        nnodes = length(mesh[!, coor[1]])
        data = zeros(dof, nnodes)
        for i in 1:dof
            data[i, :] = values(mesh[!, coor[i]])
        end

        for (name, filter) in bond_filters[2]
            if filter["Type"] == "Disk"
                nlist = disk_filter(nnodes, data, filter, nlist, dof)
            elseif filter["Type"] == "Rectangular_Plane"
                nlist = rectangular_plane_filter(nnodes, data, filter, nlist, dof)
            end
        end
        @info "Finished applying bond filters"
    end
    return nlist
end

"""
    disk_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::Vector{Vector{Int64}}, dof::Int64)

Apply the disk filter to the neighborhood list.

# Arguments
- `nnodes::Int64`: The number of nodes.
- `data::Matrix{Float64}`: The data.
- `filter::Dict`: The filter.
- `nlist::Vector{Vector{Int64}}`: The neighborhood list.
- `dof::Int64`: The degrees of freedom.
# Returns
- `nlist::Vector{Vector{Int64}}`: The filtered neighborhood list.
"""
function disk_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::Vector{Vector{Int64}}, dof::Int64)
    if dof == 2
        center = [filter["Center_X"], filter["Center_Y"]]
        normal = [filter["Normal_X"], filter["Normal_Y"]]
    else
        center = [filter["Center_X"], filter["Center_Y"], filter["Center_Z"]]
        normal = [filter["Normal_X"], filter["Normal_Y"], filter["Normal_Z"]]
    end
    #normalize vector
    normal = normal ./ norm(normal)
    for i in 1:nnodes
        filter_flag = fill(true, length(nlist[i]))
        for (jId, neighbor) in enumerate(nlist[i])
            filter_flag[jId] = !bondIntersectsDisk(data[:, i], data[:, neighbor], center, normal, filter["Radius"])
        end
        nlist[i] = nlist[i][filter_flag]
    end
    return nlist
end

"""
    rectangular_plane_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::Vector{Vector{Int64}}, dof::Int64)

Apply the rectangular plane filter to the neighborhood list.

# Arguments
- `nnodes::Int64`: The number of nodes.
- `data::Matrix{Float64}`: The data.
- `filter::Dict`: The filter.
- `nlist::Vector{Vector{Int64}}`: The neighborhood list.
- `dof::Int64`: The degrees of freedom.
# Returns
- `nlist::Vector{Vector{Int64}}`: The filtered neighborhood list.
"""
function rectangular_plane_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::Vector{Vector{Int64}}, dof::Int64)
    if dof == 2
        normal = [filter["Normal_X"], filter["Normal_Y"]]
        lower_left_corner = [filter["Lower_Left_Corner_X"], filter["Lower_Left_Corner_Y"]]
        bottom_unit_vector = [filter["Bottom_Unit_Vector_X"], filter["Bottom_Unit_Vector_Y"]]
    else
        normal = [filter["Normal_X"], filter["Normal_Y"], filter["Normal_Z"]]
        lower_left_corner = [filter["Lower_Left_Corner_X"], filter["Lower_Left_Corner_Y"], filter["Lower_Left_Corner_Z"]]
        bottom_unit_vector = [filter["Bottom_Unit_Vector_X"], filter["Bottom_Unit_Vector_Y"], filter["Bottom_Unit_Vector_Z"]]
    end
    #normalize vector
    normal = normal ./ norm(normal)
    bottom_unit_vector = bottom_unit_vector ./ norm(bottom_unit_vector)
    bottom_length = filter["Bottom_Length"]
    side_length = filter["Side_Length"]
    for iID in 1:nnodes
        filter_flag = fill(true, length(nlist[iID]))
        for (jID, neighborID) in enumerate(nlist[iID])
            intersect_inf_plane, x = bondIntersectInfinitePlane(data[:, iID], data[:, neighborID], lower_left_corner, normal)
            bond_intersect = false
            if intersect_inf_plane
                bond_intersect = bondIntersectRectanglePlane(x, lower_left_corner, bottom_unit_vector, normal, side_length, bottom_length)
            end
            filter_flag[jID] = !(intersect_inf_plane && bond_intersect)
        end
        nlist[iID] = nlist[iID][filter_flag]
    end
    return nlist
end

"""
    glob_to_loc(distribution)

Get the global to local mapping

# Arguments
- `distribution`: The distribution
# Returns
- `glob_to_loc`: The global to local mapping
"""
function glob_to_loc(distribution)
    glob_to_loc = Dict{Int64,Int64}()
    for id in eachindex(distribution)
        glob_to_loc[distribution[id]] = id
    end
    return glob_to_loc
end

end