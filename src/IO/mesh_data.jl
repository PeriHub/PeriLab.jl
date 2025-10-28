# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using LinearAlgebra
using AbaqusReader
using DataFrames
using OrderedCollections: OrderedDict
using TimerOutputs: @timeit

using ..Data_Manager
include("bond_filters.jl")
include("gcode.jl")
include("volume.jl")
using ..Helpers: fastdot, get_nearest_neighbors
using ..Logging_Module: print_table
using ..Parameter_Handling: get_mesh_name, get_header, get_node_sets,
                            get_external_topology_name, get_horizon
using ..Geometry: bond_geometry!

#export read_mesh
#export load_mesh_and_distribute
export init_data

const TOLERANCE = 1.0e-14

"""
    init_data(params::Dict, path::String, comm::MPI.Comm)

Initializes the data for the mesh.

# Arguments
- `params::Dict`: The parameters for the simulation.
- `path::String`: The path to the mesh file.
- `comm::MPI.Comm`: The MPI communicator.
# Returns
- `params::Dict`: The parameters for the simulation.
"""
function init_data(params::Dict,
                   path::String,
                   comm::MPI.Comm)
    @timeit "init_data - mesh_data,jl" begin
        size = MPI.Comm_size(comm)
        rank = MPI.Comm_rank(comm) + 1
        fem_active::Bool = false
        if rank == 1
            @timeit "load_and_evaluate_mesh" distribution,
                                             mesh,
                                             ntype,
                                             overlap_map,
                                             nlist,
                                             nlist_filtered_ids,
                                             bond_norm,
                                             dof,
                                             nsets,
                                             topology,
                                             element_distribution=load_and_evaluate_mesh(params,
                                                                                         path,
                                                                                         size,
                                                                                         Data_Manager.get_silent())
            if !isnothing(element_distribution)
                fem_active = true
            end

            data = ["Mesh input overview" "" ""
                    "Number of nodes" "."^10 length(mesh[!, "x"])
                    "Geometrical degrees of freedoms" "."^10 dof]
            if fem_active
                data = vcat(data, ["Number of finite elements" "."^10 length(topology)])
            end
            print_table(data)
        else
            dof::Int64 = 0
            distribution = nothing
            element_distribution = nothing
            mesh = DataFrame()
            nlist = nothing
            nlist_filtered_ids = nothing
            bond_norm = nothing
            nsets = nothing
            ntype = Dict("controllers" => 0, "responder" => 0)
            overlap_map = nothing
        end
        if size > 1
            @info "Synchronize data over cores"
            fem_active = broadcast_value(comm, fem_active)
            dof = broadcast_value(comm, dof)
            overlap_map = broadcast_value(comm, overlap_map)
            distribution = broadcast_value(comm, distribution)
            nsets = broadcast_value(comm, nsets)
            nlist = broadcast_value(comm, nlist)
            nlist_filtered_ids = broadcast_value(comm, nlist_filtered_ids)
        end
        num_controller::Int64 = send_single_value_from_vector(comm, 0, ntype["controllers"],
                                                              Int64)
        num_responder::Int64 = send_single_value_from_vector(comm, 0, ntype["responder"],
                                                             Int64)
        Data_Manager.set_dof(dof)
        Data_Manager.set_overlap_map(overlap_map)
        Data_Manager.set_num_controller(num_controller)
        Data_Manager.set_num_responder(num_responder)
        @debug "Get node sets"
        define_nsets(nsets)
        # defines the order of the global nodes to the local core nodes
        Data_Manager.set_distribution(distribution[rank])
        Data_Manager.set_glob_to_loc(create_global_to_local_mapping(distribution[rank]))
        @timeit "get_local_overlap_map" overlap_map=get_local_overlap_map(overlap_map,
                                                                          distribution,
                                                                          size)
        @timeit "distribution_to_cores" distribution_to_cores(comm,
                                                              mesh,
                                                              distribution,
                                                              dof)
        @timeit "distribute_neighborhoodlist_to_cores" distribute_neighborhoodlist_to_cores(comm,
                                                                                            nlist,
                                                                                            distribution,
                                                                                            false)

        if !isnothing(nlist_filtered_ids)
            create_and_distribute_bond_norm(comm,
                                            nlist_filtered_ids,
                                            distribution,
                                            bond_norm,
                                            dof)
        end

        contact_basis(params, mesh, comm, rank)

        get_bond_geometry() # gives the initial length and bond damage
        Data_Manager.set_fem(fem_active)
        if fem_active
            @debug "Set and synchronize elements"
            element_distribution = broadcast_value(comm, element_distribution)
            topology = broadcast_value(comm, topology)
            Data_Manager.set_num_elements(length(element_distribution[rank]))
            @debug "Set local topology vector"
            get_local_element_topology(topology[element_distribution[rank]],
                                       distribution[rank])
        end
        @debug "Finish init data"
    end
    barrier(comm)
    mesh = nothing
    return params
end

"""
    create_and_distribute_bond_norm(comm::MPI.Comm, nlist_filtered_ids::Vector{Vector{Int64}}, distribution::Vector{Int64}, bond_norm::Vector{Float64}, dof::Int64)

Create and distribute the bond norm

# Arguments
- `comm::MPI.Comm`: MPI communicator
- `nlist_filtered_ids::Vector{Vector{Int64}}`: The filtered neighborhood list
- `distribution::Vector{Int64}`: The distribution
- `bond_norm::Vector{Float64}`: The bond norm
- `dof::Int64`: The degree of freedom
"""
function create_and_distribute_bond_norm(comm::MPI.Comm,
                                         nlist_filtered_ids::Vector{Vector{Int64}},
                                         distribution::Vector{Vector{Int64}},
                                         bond_norm::Vector{Any},
                                         dof::Int64)
    bond_norm = broadcast_value(comm, bond_norm)
    bond_norm_field = Data_Manager.create_constant_bond_field("Bond Norm", Float64, dof, 1)
    distribute_neighborhoodlist_to_cores(comm,
                                         nlist_filtered_ids,
                                         distribution,
                                         true)
    copyto!.(bond_norm_field, bond_norm)
end

function contact_basis(params::Dict, mesh::DataFrame, comm,
                       rank::Int64)
    if !haskey(params, "Contact")
        return
    end
    ## All coordinates and block Ids are stored at all cores. Reason is, that you might need them for self contact, etc.
    if rank == 1
        mesh_id = ["x", "y"]
        if Data_Manager.get_dof() == 3
            mesh_id = ["x", "y", "z"]
        end
        points = Matrix{Float64}(mesh[:, mesh_id])
        blocks = Vector{Int64}(mesh[:, "block_id"])
    end
    if rank > 1
        blocks = nothing
        points = nothing
    end
    points = broadcast_value(comm, points)
    blocks = broadcast_value(comm, blocks)
    Data_Manager.set_all_positions(points)
    Data_Manager.set_all_blocks(blocks)
end
"""
    get_local_element_topology(topology::Vector{Vector{Int64}}, distribution::Vector{Int64})

Get the local element topology

# Arguments
- `topology::Vector{Vector{Int64}}`: The topology
- `distribution::Vector{Int64}`: The distribution
"""
function get_local_element_topology(topology::Vector{Vector{Int64}},
                                    distribution::Vector{Int64})
    if isempty(topology[1])
        return
    end
    master_len = length(topology[1][:])
    for top in topology
        if length(top) != master_len
            @error "Only one element type is supported. Please define the same numbers of nodes per element."
            return nothing
            # - new field like the bond field has to be defined for elements in the Data_Manager
            # - can be avoided right now by setting zeros in the topology vector as empty nodes
        end
    end
    topo = Data_Manager.create_constant_free_size_field("FE Topology",
                                                        Int64,
                                                        (length(topology), master_len))
    ilocal = create_global_to_local_mapping(distribution)
    for el_id in eachindex(topology)
        topo[el_id, :] = local_nodes_from_dict(ilocal, topology[el_id])
    end
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
function get_local_overlap_map(overlap_map,
                               distribution::Vector{Vector{Int64}},
                               ranks::Int64)
    if ranks == 1
        return overlap_map
    end
    for irank in 1:ranks
        ilocal = create_global_to_local_mapping(distribution[irank])
        for jrank in 1:ranks
            if irank != jrank
                overlap_map[irank][jrank]["Responder"] .= local_nodes_from_dict(ilocal,
                                                                                overlap_map[irank][jrank]["Responder"])
                overlap_map[irank][jrank]["Controller"] .= local_nodes_from_dict(ilocal,
                                                                                 overlap_map[irank][jrank]["Controller"])
            end
        end
    end
    return sort!(OrderedDict(overlap_map))
end

"""
    local_nodes_from_dict(create_global_to_local_mapping::Dict{Int,Int}, global_nodes::Vector{Int64})

Changes entries in the overlap map from the global numbering to the local computer core one.

# Arguments
- `create_global_to_local_mapping::Dict{Int,Int}`: global to local mapping
- `global_nodes::Vector{Int64}`: global nodes
# Returns
- `overlap_map::Dict{Int64, Dict{Int64, String}}`: returns overlap map with local nodes.
"""
function local_nodes_from_dict(create_global_to_local_mapping::Dict{Int,Int},
                               global_nodes::Vector{Int64})
    return Int64[create_global_to_local_mapping[global_node]
                 for
                 global_node in global_nodes
                 if haskey(create_global_to_local_mapping, global_node)]
end

"""
    distribute_neighborhoodlist_to_cores(comm::MPI.Comm, nlist, distribution)

Distributes the neighborhood list to the cores.

# Arguments
- `comm::MPI.Comm`: MPI communicator
- `nlist`: neighborhood list
- `distribution Array{Int64}`: global nodes distribution at cores
"""
function distribute_neighborhoodlist_to_cores(comm::MPI.Comm,
                                              nlist::Vector{Vector{Int64}},
                                              distribution::Vector{Vector{Int64}},
                                              filtered::Bool)
    send_msg = 0
    if filtered
        length_nlist = Data_Manager.create_constant_node_field("Number of Filtered Neighbors",
                                                               Int64, 1)
    else
        length_nlist = Data_Manager.create_constant_node_field("Number of Neighbors", Int64,
                                                               1)
    end
    rank = MPI.Comm_rank(comm)
    if rank == 0
        send_msg = get_number_of_neighbornodes(nlist, filtered)
    end
    length_nlist .= send_vector_from_root_to_core_i(comm, send_msg, length_nlist,
                                                    distribution)
    if filtered
        nlist_core = Data_Manager.create_constant_bond_field("FilteredNeighborhoodlist",
                                                             Int64, 1)
    else
        nlist_core = Data_Manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    end

    nlist_core .= nlist[distribution[rank + 1][:]]
    nlist_core .= get_local_neighbors(Data_Manager.get_local_nodes, nlist_core)
    nlist = 0
end

"""
    get_local_neighbors(mapping, nlist_core)

Gets the local neighborhood list from the global neighborhood list

# Arguments
- `mapping`: mapping function
- `nlist_core`: global neighborhood list
# Returns
- `nlist_core`: local neighborhood list
"""
function get_local_neighbors(mapping, nlist_core)
    for id in eachindex(nlist_core)
        nlist_core[id] = mapping(nlist_core[id])
    end
    return nlist_core
end

"""
    get_bond_geometry()

Gets the bond geometry
"""
function get_bond_geometry()
    dof = Data_Manager.get_dof()
    nnodes = Data_Manager.get_nnodes()
    nlist = Data_Manager.get_nlist()
    coor = Data_Manager.get_field("Coordinates")
    undeformed_bond = Data_Manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    undeformed_bond_length = Data_Manager.create_constant_bond_field("Bond Length", Float64,
                                                                     1)
    bond_geometry!(undeformed_bond,
                   undeformed_bond_length,
                   Vector{Int64}(1:nnodes),
                   nlist,
                   coor)
end

"""
    define_nsets(nsets::Dict{String,Vector{Int64}})

Defines the node sets

# Arguments
- `nsets::Dict{String,Vector{Int64}}`: Node sets read from files
"""
function define_nsets(nsets::Dict{String,Vector{Int64}})
    for nset in keys(nsets)
        Data_Manager.set_nset(nset, nsets[nset])
    end
end

"""
    distribution_to_cores(comm::MPI.Comm, mesh, distribution, dof::Int64)

Distributes the mesh data to the cores

# Arguments
- `comm::MPI.Comm`: MPI communicator
- `mesh`: Mesh
- `distribution Array{Int64}`: global nodes distribution at cores
- `dof::Int64`: Degree of freedom
"""
function distribution_to_cores(comm::MPI.Comm,
                               mesh::DataFrame,
                               distribution::Vector{Vector{Int64}},
                               dof::Int64)
    # init block_id field
    rank = MPI.Comm_rank(comm)
    if rank == 0
        meshdata = check_mesh_elements(mesh, dof)
    else
        meshdata = Dict()
    end
    meshdata = broadcast_value(comm, meshdata)
    for fieldname in keys(meshdata)
        field_dof = length(meshdata[fieldname]["Mesh ID"])
        datafield = Data_Manager.create_constant_node_field(fieldname,
                                                            meshdata[fieldname]["Type"],
                                                            field_dof)
        for (localDof, mesh_id) in enumerate(meshdata[fieldname]["Mesh ID"])
            if rank == 0
                send_msg = meshdata[fieldname]["Type"].(mesh[!, mesh_id])
                # example send_msg = Float64.(mesh[!, names(mesh)[idof]])
                # redefine from Float64 standard to Float64 for MPI
            else
                send_msg = 0
            end
            if field_dof == 1
                datafield .= send_vector_from_root_to_core_i(comm, send_msg, datafield,
                                                             distribution)
            else
                datafield[:,
                localDof] = send_vector_from_root_to_core_i(comm,
                                                                         send_msg,
                                                                         datafield[:,
                                                                         localDof],
                                                                         distribution)
            end
        end
    end
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
    mesh_info_dict = Dict{String,Dict{String,Any}}()

    for (id, mesh_entry) in enumerate(mnames)
        field_dof = 1
        if ("y" == mesh_entry) ||
           ("z" == mesh_entry) ||
           (mesh_entry[1:(end - 1)] in keys(mesh_info_dict))
            continue
        end

        if "x" == mesh_entry
            name = "Coordinates"
            mesh_id = ["x", "y"]
            if dof == 3
                mesh_id = ["x", "y", "z"]
            end
        elseif "volume" == mesh_entry
            name = "Volume"
            mesh_id = ["volume"]
        elseif "block_id" == mesh_entry
            name = "Block_Id"
            mesh_id = ["block_id"]
        else
            if "_x" == mesh_entry[(end - 1):end]
                if id + 1 <= length(mnames)
                    if mnames[id + 1][(end - 1):end] == "_y"
                        name = mesh_entry[1:(end - 2)]
                        mesh_id = [name * "_x", name * "_y"]
                        if dof == 3
                            mesh_id = [name * "_x", name * "_y", name * "_z"]
                        end
                    end
                end
            elseif mesh_entry[(end - 1):end] in ["_y", "_z"]
                continue
            else
                name = mesh_entry
                mesh_id = [name]
            end
        end

        if mesh[1, mesh_id[1]] isa Bool
            vartype = Bool
        elseif name in ["Coordinates", "Volume"]
            vartype = Float64
        else
            vartype = typeof(sum(sum(mesh[:, mid]) for mid in mesh_id))
        end
        mesh_info_dict[name] = Dict{String,Any}("Mesh ID" => mesh_id, "Type" => vartype)
    end
    if !(haskey(mesh_info_dict, "Coordinates"))
        @error "No coordinates defined"
        return nothing
    end
    if !(haskey(mesh_info_dict, "Block_Id"))
        @error "No blocks defined"
        return nothing
    end
    if !(haskey(mesh_info_dict, "Volume"))
        @error "No volumes defined"
        return nothing
    end
    return mesh_info_dict
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
    return csv_reader(filename)
end

"""
    csv_reader(filename::String)

Read csv and return it as a DataFrame.

# Arguments
- `filename::String`: The path to the mesh file.
# Returns
- `csvData::DataFrame`: The csv data a DataFrame.
"""
function csv_reader(filename::String)
    header_line, header = get_header(filename)
    return CSV.read(filename,
                    DataFrame;
                    delim = " ",
                    ignorerepeated = true,
                    header = header,
                    skipto = header_line + 1,
                    comment = "#",)
end

"""
    read_mesh(filename::String, params::Dict)

Read mesh data from a file and return it as a DataFrame.

# Arguments
- `filename::String`: The path to the mesh file.
- `params::Dict`: The input parameters.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
"""
function read_mesh(filename::String, params::Dict)
    if !isfile(filename)
        @error "File $filename does not exist"
        return nothing
    end

    @info "Read mesh file $filename"

    if params["Discretization"]["Type"] == "Exodus"
        exo = ExodusDatabase(filename, "r")

        coords = read_coordinates(exo)
        mesh_df = DataFrame(x = Float64[],
                            y = Float64[],
                            z = Float64[],
                            volume = Float64[],
                            block_id = Int64[])
        block_ids = read_ids(exo, Block)

        for (iID, block_id) in enumerate(block_ids)
            block = read_block(exo, block_id)
            block_id_map = Exodus.read_block_connectivity(exo,
                                                          block_id,
                                                          block.num_nodes_per_elem *
                                                          block.num_elem)
            if block.elem_type == "TETRA"
                for i in 1:(block.num_elem)
                    indices = (block.num_nodes_per_elem * (i - 1) + 1):(block.num_nodes_per_elem * i)
                    node_ids = block_id_map[indices]
                    vertices = coords[:, node_ids]
                    center = sum(vertices, dims = 2) / size(vertices)[2]
                    volume = tetrahedron_volume(vertices)
                    push!(mesh_df,
                          (x = center[1],
                           y = center[2],
                           z = center[3],
                           volume = volume,
                           block_id = Int64(block_id)))
                end
            elseif block.elem_type == "HEX8"
                for i in 1:(block.num_elem)
                    indices = (block.num_nodes_per_elem * (i - 1) + 1):(block.num_nodes_per_elem * i)
                    node_ids = block_id_map[indices]
                    vertices = coords[:, node_ids]
                    center = sum(vertices, dims = 2) / size(vertices)[2]
                    volume = hex8_volume(vertices)
                    push!(mesh_df,
                          (x = center[1],
                           y = center[2],
                           z = center[3],
                           volume = volume,
                           block_id = Int64(block_id)))
                end
            else
                @error "Element type $(block.elem_type) not supported"
            end
        end

        close(exo)
        coords = nothing
        block_ids = nothing

        return mesh_df

    elseif params["Discretization"]["Type"] == "Abaqus"
        mesh = abaqus_read_mesh(filename; verbose = false)

        nodes = mesh["nodes"]
        elements = mesh["elements"]
        element_sets = mesh["element_sets"]
        element_types = mesh["element_types"]

        dof = 2
        nodes_vector = collect(values(nodes))
        if size(nodes_vector[1])[1] == 3
            dof = 3
        end
        @info "Abaqus mesh with $dof DOF"

        num_elements = length(elements)
        mesh_df = ifelse(dof == 2,
                         DataFrame(x = Array{Float64,1}(undef, num_elements),
                                   y = Array{Float64,1}(undef, num_elements),
                                   volume = Array{Float64,1}(undef, num_elements),
                                   block_id = Array{Int64,1}(undef, num_elements)),
                         DataFrame(x = Array{Float64,1}(undef, num_elements),
                                   y = Array{Float64,1}(undef, num_elements),
                                   z = Array{Float64,1}(undef, num_elements),
                                   volume = Array{Float64,1}(undef, num_elements),
                                   block_id = Array{Int64,1}(undef, num_elements)))

        id = 1
        block_id = 1
        # element_written = Array{Int64,1}(undef, num_elements)
        element_written = []
        nsets = Dict{String,Vector{Int64}}()

        nset_names = []

        for boundary_condtion in keys(params["Boundary Conditions"])
            if haskey(params["Boundary Conditions"][boundary_condtion], "Node Set")
                push!(nset_names,
                      params["Boundary Conditions"][boundary_condtion]["Node Set"])
            end
        end
        nset_names = unique(nset_names)

        # sort element_sets by length
        # element_sets_keys = sort(collect(keys(element_sets)), by=x -> length(element_sets[x]), rev=true)
        element_sets_keys = collect(keys(element_sets))
        for nset in nset_names
            if nset in element_sets_keys
                deleteat!(element_sets_keys, findfirst(x -> x == nset, element_sets_keys))
                push!(element_sets_keys, nset)
            end
        end
        block_names = []
        for key in element_sets_keys
            element_set = element_sets[key]
            ns_nodes = Array{Int64,1}(undef, length(element_set))
            nset_only = true
            for (jID, element_id) in enumerate(element_set)
                if element_id in element_written
                    # push!(ns_nodes, findfirst(x -> x == element_id, element_written))
                    if key in nset_names
                        ns_nodes[jID] = findfirst(x -> x == element_id, element_written)
                    end
                    continue
                end
                nset_only = false
                ns_nodes[jID] = id
                node_ids = elements[element_id]
                element_type = element_types[element_id]
                vertices = [nodes[node_id] for node_id in node_ids]
                volume = calculate_volume(string(element_type), vertices)
                center = sum(vertices) / length(vertices)
                if dof == 2
                    mesh_df[id, :] = [center[1], center[2], volume, block_id]
                else
                    mesh_df[id, :] = [center[1], center[2], center[3], volume, block_id]
                end
                # element_written[id] = element_id
                push!(element_written, element_id)
                id += 1
            end
            if key in nset_names
                nsets[key] = ns_nodes
            end
            if !nset_only
                block_id += 1
                push!(block_names, key)
            end
        end
        @info "Found $(block_id-1) block(s)"
        @info "Blocks: $block_names"
        @info "Found $(length(nsets)) node set(s)"
        @info "NodeSets: $(keys(nsets))"

        mesh = nothing
        nodes = nothing
        elements = nothing
        element_sets = nothing

        return mesh_df, nsets

    elseif params["Discretization"]["Type"] in ["Text File", "Gcode"]
        return csv_reader(filename)
    else
        @error "Discretization type not supported"
        return nothing
    end
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
    @info "2d problem assumed, please define plane stress or plane strain if needed in Models"
    return Int64(2)
end

"""
    check_for_duplicate_in_dataframe(mesh::DataFrame)

check duplicated entries and throws an error if one is there. If not everything is ok.

# Arguments
- `mesh::DataFrame`: The input mesh data represented as a DataFrame.
"""
function check_for_duplicate_in_dataframe(mesh::DataFrame)
    duplicates = findall(nonunique(mesh))
    if length(duplicates) > 0
        @error "Mesh contains duplicate nodes! Nodes: $duplicates"
        return true
    end
    return false
end

"""
    check_types_in_dataframe(mesh::DataFrame)

check if block_id in mesh contains only int.

# Arguments
- `mesh::DataFrame`: The input mesh data represented as a DataFrame.
"""
function check_types_in_dataframe(mesh::DataFrame)
    # check if block_id in mesh contains only int
    if !(eltype(mesh.block_id) <: Integer)
        @error "block_id in mesh is $(eltype(mesh.block_id)), but it should be an Integer!"
        return false
    end
    return true
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
- `nsets::Dict`: The node sets
- `topology::Int64`::Array{Int64,nelement:nodes}`: The topology of elements.
- `el_distribution::Array{Int64,1}`: The distribution of the finite elements.
"""
function load_and_evaluate_mesh(params::Dict,
                                path::String,
                                ranksize::Int64,
                                silent::Bool)
    filename = joinpath(path, get_mesh_name(params))
    if params["Discretization"]["Type"] == "Abaqus"
        @timeit "read_mesh" mesh, nsets=read_mesh(filename, params)
    elseif params["Discretization"]["Type"] == "Gcode"
        txt_file = replace(filename, ".gcode" => ".txt")
        @info txt_file
        if params["Discretization"]["Gcode"]["Overwrite Mesh"] || !isfile(txt_file)
            mesh = get_gcode_mesh(filename, params, silent)
        else
            mesh = read_mesh(txt_file, params)
        end
        nsets = get_node_sets(params, path, mesh)
    else
        @debug "Read node sets"
        @timeit "read_mesh" mesh=read_mesh(filename, params)
        nsets = get_node_sets(params, path, mesh)
    end
    nnodes = size(mesh, 1) + 1
    mesh, surface_ns = extrude_surface_mesh(mesh, params)
    if !isnothing(surface_ns)
        for (key, values) in surface_ns
            nsets[key] = Vector{Int64}(values .+ nnodes)
        end
    end
    check_for_duplicate_in_dataframe(mesh)
    check_types_in_dataframe(mesh)

    external_topology = nothing
    if !isnothing(get_external_topology_name(params, path))
        external_topology = read_external_topology(joinpath(path,
                                                            get_external_topology_name(params,
                                                                                       path)))
    end
    if !isnothing(external_topology)
        @info "External topology files was read."
    end
    dof::Int64 = set_dof(mesh)
    @timeit "neighborhoodlist" nlist, _=create_neighborhoodlist(mesh, params, dof)
    @debug "Finished init Neighborhoodlist"
    @timeit "apply_bond_filters" nlist, nlist_filtered_ids,
                                 bond_norm=apply_bond_filters(nlist,
                                                              mesh,
                                                              params,
                                                              dof)
    topology = nothing
    if !isnothing(external_topology)
        @info "Create a consistent neighborhood list with external topology definition."
        nlist,
        topology = create_consistent_neighborhoodlist(external_topology,
                                                      params["Discretization"]["Input External Topology"],
                                                      nlist,
                                                      dof)
    end
    @debug "Start distribution"
    if haskey(params["Discretization"], "Distribution Type")
        @timeit "node_distribution" distribution, ptc,
                                    ntype=node_distribution(nlist,
                                                            ranksize,
                                                            params["Discretization"]["Distribution Type"])
    else
        @timeit "node_distribution" distribution, ptc,
                                    ntype=node_distribution(nlist,
                                                            ranksize)
    end

    el_distribution = nothing
    if haskey(params, "FEM") && !isnothing(external_topology)
        @debug "Start element distribution"
        el_distribution = element_distribution(topology, ptc, ranksize)
    end
    @debug "Finished distribution"
    @debug "Create Overlap"
    @timeit "overlap_map" overlap_map=create_overlap_map(distribution, ptc, ranksize)
    @debug "Finished Overlap"
    return distribution,
           mesh,
           ntype,
           overlap_map,
           nlist,
           nlist_filtered_ids,
           bond_norm,
           dof,
           nsets,
           topology,
           el_distribution
end

function create_consistent_neighborhoodlist(external_topology::DataFrame,
                                            params::Dict,
                                            nlist::Vector{Vector{Int64}},
                                            dof::Int64)
    pd_neighbors::Bool = false
    if haskey(params, "Add Neighbor Search")
        pd_neighbors = params["Add Neighbor Search"]
    end
    number_of_elements = length(external_topology[:, 1])
    topology::Vector{Vector{Int64}} = []
    for i_el in 1:number_of_elements
        push!(topology, collect(skipmissing(external_topology[i_el, :])))
    end
    # creates a field of length highest point number occuring in finite element
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
        if !pd_neighbors # deletes the FE neighbor points from the PD points
            for nID in nlist[fe_node]
                index = findfirst(x -> x == fe_node, nlist[nID])
                if !isnothing(index)
                    deleteat!(nlist[nID], index)
                end
            end
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
    return neighbors(mesh, params, coor[1:dof])
end

"""
    get_number_of_neighbornodes(nlist::Vector{Vector{Int64}})

Get the number of neighbors for each node.

# Arguments
- `nlist::Vector{Vector{Int64}}`: The neighborhood list of the mesh elements.
# Returns
- `length_nlist::Vector{Int64}`: The number of neighbors for each node.
"""
function get_number_of_neighbornodes(nlist::Vector{Vector{Int64}}, filtered::Bool)
    len = length(nlist)
    length_nlist = zeros(Int64, len)
    for id in 1:len
        if !filtered && length(nlist[id]) == 0
            @error "Node $id has no neighbors please check the horizon."
            return nothing
        end
        length_nlist[id] = length(nlist[id])
    end
    return length_nlist
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
function element_distribution(topology::Vector{Vector{Int64}},
                              ptc::Vector{Int64},
                              size::Int64)
    nelements = length(topology)
    if size == 1
        distribution = [collect(1:nelements)]
        etc::Vector{Int64} = []
    else
        distribution, etc = create_distribution(nelements, size)
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
                    distribution[etc[el_id]] = filter(x -> x != el_id,
                                                      distribution[etc[el_id]])
                    etc[el_id] = min_core
                end
            end
        end
    end
    return distribution
end

"""
    node_distribution(nlist::Vector{Vector{Int64}}, size::Int64)

Create the distribution of the nodes.

# Arguments
- `nlist::Vector{Vector{Int64}}`: The neighborhood list of the mesh elements.
- `size::Int64`: The number of ranks.
- `distribution_type::String`: The distribution type.
# Returns
- `distribution::Vector{Vector{Int64}}`: The distribution of the nodes.
- `ptc::Vector{Int64}`: Defines at which core / rank each node lies.
- `ntype::Dict`: The type of the nodes.
"""
function node_distribution(nlist::Vector{Vector{Int64}},
                           size::Int64,
                           distribution_type::String = "Neighbor based")
    nnodes = length(nlist)
    ntype = Dict("controllers" => Int64[], "responder" => Int64[])
    if size == 1
        append!(ntype["controllers"], nnodes)
        append!(ntype["responder"], 0)
        return [collect(1:nnodes)], Int64[], ntype
    end

    if distribution_type == "Neighbor based"
        distribution, ptc = create_distribution_neighbor_based(nnodes, nlist, size)
    elseif distribution_type == "Node based"
        distribution, ptc = create_distribution_node_based(nnodes, nlist, size)
    else
        distribution, ptc = create_distribution(nnodes, size)
    end
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
        append!(distribution[i], sort!(unique(tempid)))
        append!(ntype["responder"], length(unique(tempid)))
    end
    return distribution, ptc, ntype
end

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
                overlap_map[i][j] = Dict{String,Vector{Int64}}("Responder" => Int64[],
                                                               "Controller" => Int64[])
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
function create_overlap_map(distribution::Vector{Vector{Int64}},
                            ptc::Vector{Int64},
                            size::Int64)
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
    create_distribution_node_based(nnodes::Int64,nlist::Vector{Vector{Int64}}, size::Int64)

Calculate the initial size of each chunk for a nearly equal number of nodes vs. cores this algorithm might lead to the problem, that the last core is not equally loaded

# Arguments
- `nnodes::Int64`: The number of nodes.
- `nlist::Vector{Vector{Int64}}`: The neighborhood list.
- `size::Int64`: The number of cores.
# Returns
- `distribution::Array{Int64,1}`: The distribution of the nodes.
- `point_to_core::Array{Int64,1}`: The number of nodes in each rank.
"""
function create_distribution_node_based(nnodes::Int64,
                                        nlist::Vector{Vector{Int64}},
                                        size::Int64)
    if size > nnodes
        @error "Number of cores $size exceeds number of nodes $nnodes."
        return nothing, nothing
    end
    chunk_size = div(nnodes, size)
    if size == 1
        return create_distribution(nnodes, size)
    end
    # Split the data into chunks
    distribution = Vector{Vector{Int64}}(undef, size - 1)
    for i in 1:(size - 1)
        distribution[i] = zeros(Int64, chunk_size)
    end
    # this leads to a coupling of elements
    #distribution::Vector{Vector{Int64}} = fill(Vector{Int64}(zeros(Int64, chunk_size)), size - 1)
    push!(distribution, zeros(Int64, chunk_size + nnodes - size * chunk_size))

    point_to_core::Vector{Int64} = zeros(Int64, nnodes)
    not_included_nodes::Vector{Bool} = fill(true, nnodes)
    isize::Int64 = 1

    for iID in 1:nnodes
        if not_included_nodes[iID]
            idx = findfirst(isequal(0), distribution[isize])
            if isnothing(idx)
                isize += 1
                idx = 1
            end
            distribution[isize][idx] = iID
            not_included_nodes[iID] = false
        end # muss das hinter die nächste schleife? -> test über die bool liste
        for jID in nlist[iID]
            if not_included_nodes[jID]
                idx = findfirst(isequal(0), distribution[isize])
                if isnothing(idx)
                    isize += 1
                    idx = 1
                end
                distribution[isize][idx] = jID
                not_included_nodes[jID] = false
            end
        end
    end
    if sum(not_included_nodes) > 0
        @error "code must be improved"
    end
    for i in 1:size
        point_to_core[distribution[i]] .= i
    end
    return distribution, point_to_core
end

"""
    create_distribution_neighbor_based(nnodes::Int64,nlist::Vector{Vector{Int64}}, size::Int64)

Calculate the initial size of each chunk for a nearly equal number of nodes vs. cores this algorithm might lead to the problem, that the last core is not equally loaded

# Arguments
- `nnodes::Int64`: The number of nodes.
- `nlist::Vector{Vector{Int64}}`: The neighborhood list.
- `size::Int64`: The number of cores.
# Returns
- `distribution::Array{Int64,1}`: The distribution of the nodes.
- `point_to_core::Array{Int64,1}`: The number of nodes in each rank.
"""
function create_distribution_neighbor_based(nnodes::Int64,
                                            nlist::Vector{Vector{Int64}},
                                            size::Int64)
    if size > nnodes
        @error "Number of cores $size exceeds number of nodes $nnodes."
        return nothing, nothing
    end

    number_neighbors_sum = sum(length.(nlist))
    number_neighbors_per_core = div(number_neighbors_sum, size)

    if size == 1
        return create_distribution(nnodes, size)
    end
    # chunk_size = div(nnodes, size)
    # Split the data into chunks
    distribution = Vector{Vector{Int64}}()
    distr_sub_vector = Vector{Int64}()
    # for i in 1:size-1
    #     distribution[i] = zeros(Int64, chunk_size)
    # end
    point_to_core::Vector{Int64} = zeros(Int64, nnodes)
    not_included_nodes::Vector{Bool} = fill(true, nnodes)
    number_neighbors::Int64 = 0

    for iID in 1:nnodes
        if not_included_nodes[iID]
            number_neighbors += length(nlist[iID])
            push!(distr_sub_vector, iID)
            not_included_nodes[iID] = false
            if number_neighbors > number_neighbors_per_core
                push!(distribution, distr_sub_vector)
                distr_sub_vector = []
                number_neighbors = 0
            end
        end # muss das hinter die nächste schleife? -> test über die bool liste
        for jID in nlist[iID]
            if not_included_nodes[jID]
                number_neighbors += length(nlist[jID])
                push!(distr_sub_vector, jID)
                not_included_nodes[jID] = false
                if number_neighbors > number_neighbors_per_core
                    push!(distribution, distr_sub_vector)
                    distr_sub_vector = []
                    number_neighbors = 0
                end
            end
        end
    end
    push!(distribution, distr_sub_vector)
    if sum(not_included_nodes) > 0 || length(distribution) != size
        @error "code must be improved"
    end
    for i in 1:size
        point_to_core[distribution[i]] .= i
    end
    return distribution, point_to_core
end

"""
    create_distribution(nnodes::Int64, size::Int64)

Calculate the initial size of each chunk for a nearly equal number of nodes vs. cores this algorithm might lead to the problem, that the last core is not equally loaded

# Arguments
- `nnodes::Int64`: The number of nodes.
- `size::Int64`: The number of cores.
# Returns
- `distribution::Array{Int64,1}`: The distribution of the nodes.
- `point_to_core::Array{Int64,1}`: The number of nodes in each rank.
"""
function create_distribution(nnodes::Int64, size::Int64)
    if size > nnodes
        @error "Number of cores $size exceeds number of nodes $nnodes."
        return nothing, nothing
    end
    chunk_size = div(nnodes, size)
    # Split the data into chunks
    distribution = fill(Int64[], size)
    point_to_core = zeros(Int64, nnodes)
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
function neighbors(mesh::DataFrame, params::Dict,
                   coor::Union{Vector{Int64},Vector{String}})
    @info "Init Neighborhoodlist"
    nnodes = length(mesh[!, coor[1]])
    dof = length(coor)
    data = zeros(nnodes, dof)
    neighborList = fill(Vector{Int64}([]), nnodes)

    for i in 1:dof
        data[:, i] = values(mesh[!, coor[i]])
    end
    block_ids = unique(mesh[!, "block_id"])
    # TODO include mesh horizon
    radius = zeros(Float64, nnodes)
    for iID in 1:nnodes
        radius[iID] = get_horizon(params, mesh[!, "block_id"][iID])
    end

    return get_nearest_neighbors(1:nnodes, dof, data, data, radius, neighborList)
end

"""
    glob_to_loc(distribution)

Get the global to local mapping

# Arguments
- `distribution`: The distribution
# Returns
- `create_global_to_local_mapping`: The global to local mapping
"""
function create_global_to_local_mapping(distribution)
    create_global_to_local_mapping = Dict{Int64,Int64}()
    for id in eachindex(distribution)
        create_global_to_local_mapping[distribution[id]] = id
    end
    return create_global_to_local_mapping
end

"""
    extrude_surface_mesh(mesh::DataFrame)

extrude the mesh at the surface of the block

# Arguments
- `mesh::DataFrame`: The input mesh data represented as a DataFrame.
- `params::Dict`: The input parameters.
"""
function extrude_surface_mesh(mesh::DataFrame, params::Dict)
    if !("Surface Extrusion" in keys(params["Discretization"]))
        return mesh, nothing
    end
    direction = params["Discretization"]["Surface Extrusion"]["Direction"]
    step_x = params["Discretization"]["Surface Extrusion"]["Step_X"]
    step_y = params["Discretization"]["Surface Extrusion"]["Step_Y"]
    step_z = params["Discretization"]["Surface Extrusion"]["Step_Z"]
    number = params["Discretization"]["Surface Extrusion"]["Number"]

    # Finding min and max values for each dimension
    min_x, max_x = extrema(mesh.x)
    min_y, max_y = extrema(mesh.y)
    min_z = 0.0
    max_z = 0.0
    dof = 2
    if "z" in names(mesh)
        min_z, max_z = extrema(mesh.z)
        dof = 3
    end

    if direction == "X"
        coord_min = min_x
        coord_max = max_x
        row_min = min_y
        row_max = max_y
    elseif direction == "Y"
        coord_min = min_y
        coord_max = max_y
        row_min = min_x
        row_max = max_x
    end

    block_id = maximum(mesh.block_id) + 1

    id = 0

    node_sets = Dict("Extruded_1" => [], "Extruded_2" => [])

    for i in (coord_max + step_x):step_x:(coord_max + step_x * number),
        j in row_min:step_y:(row_max + step_y),
        k in min_z:step_z:max_z
        if direction == "X"
            if dof == 2
                push!(mesh, (x = i, y = j, volume = step_x * step_y, block_id = block_id))
            else
                push!(mesh,
                      (x = i,
                       y = j,
                       z = k,
                       volume = step_x * step_y * step_z,
                       block_id = block_id))
            end
        elseif direction == "Y"
            if dof == 2
                push!(mesh, (x = j, y = i, volume = step_x * step_y, block_id = block_id))
            else
                push!(mesh,
                      (x = j,
                       y = i,
                       z = k,
                       volume = step_x * step_y * step_z,
                       block_id = block_id))
            end
        end
        append!(node_sets["Extruded_1"], [Int64(id)])
        id += 1
    end

    block_id += 1

    for i in (coord_min - step_x):(-step_x):(coord_min - step_x * number),
        j in row_min:step_y:(row_max + step_y),
        k in min_z:step_z:max_z
        if direction == "X"
            if dof == 2
                push!(mesh, (x = i, y = j, volume = step_x * step_y, block_id = block_id))
            else
                push!(mesh,
                      (x = i,
                       y = j,
                       z = k,
                       volume = step_x * step_y * step_z,
                       block_id = block_id))
            end
        elseif direction == "Y"
            if dof == 2
                push!(mesh, (x = j, y = i, volume = step_x * step_y, block_id = block_id))
            else
                push!(mesh,
                      (x = j,
                       y = i,
                       z = k,
                       volume = step_x * step_y * step_z,
                       block_id = block_id))
            end
        end
        append!(node_sets["Extruded_2"], [Int64(id)])
        id += 1
    end
    return mesh, node_sets
end
