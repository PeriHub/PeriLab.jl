# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using LinearAlgebra
using AbaqusReader
using DataFrames
using Delaunay
using NearestNeighbors: BallTree
using NearestNeighbors: inrange
using PrettyTables
include("./logging.jl")
using .Logging_module: print_table

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
        size = MPI.Comm_size(comm)
        rank = MPI.Comm_rank(comm) + 1
        fem_active::Bool = false
        if rank == 1
            @timeit to "load_and_evaluate_mesh" distribution, mesh, ntype, overlap_map, nlist, nlist_filtered_ids, bond_norm, dof, nsets, topology, element_distribution = load_and_evaluate_mesh(params, path, size, to)
            if !isnothing(element_distribution)
                fem_active = true
            end

            data = [
                "Mesh input overview" "" "";
                "Number of nodes" "."^10 length(mesh[!, "x"]);
                "Geometrical degrees of freedoms" "."^10 dof
            ]
            if fem_active
                data = vcat(data, ["Number of finite elements" "."^10 length(topology)])
            end
            print_table(data, datamanager)
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
            fem_active = send_value(comm, 0, fem_active)
            dof = send_value(comm, 0, dof)
            overlap_map = send_value(comm, 0, overlap_map)
            distribution = send_value(comm, 0, distribution)
            nsets = send_value(comm, 0, nsets)
            nlist = send_value(comm, 0, nlist)
            nlist_filtered_ids = send_value(comm, 0, nlist_filtered_ids)
        end
        num_controller::Int64 = send_single_value_from_vector(comm, 0, ntype["controllers"], Int64)
        num_responder::Int64 = send_single_value_from_vector(comm, 0, ntype["responder"], Int64)
        dof = datamanager.set_dof(dof)
        datamanager.set_overlap_map(overlap_map)
        datamanager.set_num_controller(num_controller)
        datamanager.set_num_responder(num_responder)
        @debug "Get node sets"
        datamanager = define_nsets(nsets, datamanager)
        # defines the order of the global nodes to the local core nodes
        datamanager.set_distribution(distribution[rank])
        datamanager.set_glob_to_loc(glob_to_loc(distribution[rank]))
        @timeit to "get_local_overlap_map" overlap_map = get_local_overlap_map(overlap_map, distribution, size)
        @timeit to "distribution_to_cores" datamanager = distribution_to_cores(comm, datamanager, mesh, distribution, dof)
        @timeit to "distribute_neighborhoodlist_to_cores" datamanager = distribute_neighborhoodlist_to_cores(comm, datamanager, nlist, distribution, false)

        if !isnothing(nlist_filtered_ids)
            create_and_distribute_bond_norm(comm, datamanager, nlist_filtered_ids, distribution, bond_norm, dof)
        end

        datamanager.set_block_list(datamanager.get_field("Block_Id"))
        datamanager = get_bond_geometry(datamanager) # gives the initial length and bond damage
        datamanager.set_fem(fem_active)
        if fem_active
            @debug "Set and synchronize elements"
            element_distribution = send_value(comm, 0, element_distribution)
            topology = send_value(comm, 0, topology)
            datamanager.set_num_elements(length(element_distribution[rank]))
            @debug "Set local topology vector"
            datamanager = get_local_element_topology(datamanager, topology[element_distribution[rank]], distribution[rank])
        end
        mesh = nothing
        @debug "Finish init data"
    end
    return datamanager, params
end

"""
    create_and_distribute_bond_norm(comm::MPI.Comm, datamanager::Module, nlist_filtered_ids::Vector{Vector{Int64}}, distribution::Vector{Int64}, bond_norm::Vector{Float64}, dof::Int64)

Create and distribute the bond norm

# Arguments
- `comm::MPI.Comm`: MPI communicator
- `datamanager::Module`: Data manager
- `nlist_filtered_ids::Vector{Vector{Int64}}`: The filtered neighborhood list
- `distribution::Vector{Int64}`: The distribution
- `bond_norm::Vector{Float64}`: The bond norm
- `dof::Int64`: The degree of freedom
"""
function create_and_distribute_bond_norm(comm::MPI.Comm, datamanager::Module, nlist_filtered_ids::Vector{Vector{Int64}}, distribution::Vector{Vector{Int64}}, bond_norm::Vector{Matrix{Float64}}, dof::Int64)
    bond_norm = send_value(comm, 0, bond_norm)
    bond_norm_field = datamanager.create_constant_bond_field("Bond Norm", Float64, dof, 1)
    datamanager = distribute_neighborhoodlist_to_cores(comm, datamanager, nlist_filtered_ids, distribution, true)
    bond_norm_field .= bond_norm
end

"""
    get_local_element_topology(datamanager::Module, topology::Vector{Vector{Int64}}, distribution::Vector{Int64})

Get the local element topology

# Arguments
- `datamanager::Module`: The datamanager
- `topology::Vector{Vector{Int64}}`: The topology
- `distribution::Vector{Int64}`: The distribution
# Returns
- `datamanager::Module`: The datamanager
"""
function get_local_element_topology(datamanager::Module, topology::Vector{Vector{Int64}}, distribution::Vector{Int64})
    if length(topology[1]) == 0
        return datamanager
    end
    master_len = length(topology[1][:])
    for top in topology
        if length(top) != master_len
            @error "Only one element type is supported. Please define the same numbers of nodes per element."
            return nothing
            """
            - new field like the bond field has to be defined for elements in the datamanager
            - can be avoided right now by setting zeros in the topology vector as empty nodes
            """
        end
    end
    topo = datamanager.create_constant_free_size_field("FE Topology", Int64, (length(topology), master_len))
    ilocal = glob_to_loc(distribution)
    for el_id in eachindex(topology)
        topo[el_id, :] = local_nodes_from_dict(ilocal, topology[el_id])
    end
    return datamanager
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
                overlap_map[irank][jrank]["Responder"] .= local_nodes_from_dict(ilocal, overlap_map[irank][jrank]["Responder"])
                overlap_map[irank][jrank]["Controller"] .= local_nodes_from_dict(ilocal, overlap_map[irank][jrank]["Controller"])
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
function distribute_neighborhoodlist_to_cores(comm::MPI.Comm, datamanager::Module, nlist::Vector{Vector{Int64}}, distribution::Vector{Vector{Int64}}, filtered::Bool)
    send_msg = 0
    if filtered
        lenNlist = datamanager.create_constant_node_field("Number of Filtered Neighbors", Int64, 1)
    else
        lenNlist = datamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    end
    rank = MPI.Comm_rank(comm)
    if rank == 0
        send_msg = get_number_of_neighbornodes(nlist, filtered)
    end
    lenNlist .= send_vector_from_root_to_core_i(comm, send_msg, lenNlist, distribution)
    if filtered
        nlistCore = datamanager.create_constant_bond_field("FilteredNeighborhoodlist", Int64, 1)
    else
        nlistCore = datamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    end

    nlistCore .= nlist[distribution[rank+1][:]]
    nlistCore .= get_local_neighbors(datamanager.get_local_nodes, nlistCore)
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
    undeformed_bond = datamanager.create_constant_bond_field("Bond Geometry", Float64, dof)
    undeformed_bond_length = datamanager.create_constant_bond_field("Bond Length", Float64, 1)
    undeformed_bond, undeformed_bond_length = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, undeformed_bond, undeformed_bond_length)
    return datamanager
end

"""
    define_nsets(nsets::Dict{String,Vector{Any,Int64}}, datamanager::Module)

Defines the node sets

# Arguments
- `nsets::Dict{String,Vector{Any,Int64}}`: Node sets read from files
- `datamanager::Module`: Data manager
"""
function define_nsets(nsets::Dict{String,Any}, datamanager::Module)
    for nset in keys(nsets)
        datamanager.set_nset(nset, nsets[nset])
    end
    return datamanager
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
                datafield .= send_vector_from_root_to_core_i(comm, send_msg, datafield, distribution)
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
        return nothing
    end
    if !(haskey(meshInfoDict, "Block_Id"))
        @error "No blocks defined"
        return nothing
    end
    if !(haskey(meshInfoDict, "Volume"))
        @error "No volumes defined"
        return nothing
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
    return CSV.read(filename, DataFrame; delim=" ", ignorerepeated=true, header=header, skipto=header_line + 1, comment="#")
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
        mesh_df = DataFrame(
            x=[],
            y=[],
            z=[],
            volume=[],
            block_id=[],
        )
        block_ids = read_ids(exo, Block)

        for (iID, block_id) in enumerate(block_ids)
            block = read_block(exo, block_id)
            block_id_map = Exodus.read_block_connectivity(exo, block_id, block.num_nodes_per_elem * block.num_elem)
            if block.elem_type == "TETRA"
                for i in 1:block.num_elem
                    indices = block.num_nodes_per_elem*(i-1)+1:block.num_nodes_per_elem*i
                    node_ids = block_id_map[indices]
                    vertices = coords[:, node_ids]
                    center = (vertices[:, 1] + vertices[:, 2] + vertices[:, 3] + vertices[:, 4]) / 4.0
                    volume = abs(dot(vertices[:, 1] - vertices[:, 4], cross(vertices[:, 2] - vertices[:, 4], vertices[:, 3] - vertices[:, 4]))) / 6.0
                    push!(mesh_df, (x=center[1], y=center[2], z=center[3], volume=volume, block_id=Int64(block_id)))
                end
            end
        end

        close(exo)
        coords = nothing
        block_ids = nothing

        return mesh_df

    elseif params["Discretization"]["Type"] == "Abaqus"

        mesh = abaqus_read_mesh(filename; verbose=false)

        nodes = mesh["nodes"]
        elements = mesh["elements"]
        element_sets = mesh["element_sets"]
        element_types = mesh["element_types"]

        dof = 2

        nodes_vector = collect(values(nodes))
        node_z = [node[end] for node in nodes_vector]
        if length(unique(node_z)) > 2
            dof = 3
        end
        @info "Abaqus mesh with $dof DOF"

        mesh_df = ifelse(dof == 2,
            DataFrame(x=[], y=[], volume=[], block_id=Int[]),
            DataFrame(x=[], y=[], z=[], volume=[], block_id=Int[])
        )

        id = 1
        block_id = 1
        element_written = []
        nsets = Dict{String,Any}()

        # sort element_sets by length
        # element_sets_keys = sort(collect(keys(element_sets)), by=x -> length(element_sets[x]), rev=true)
        element_sets_keys = collect(keys(element_sets))
        for key in element_sets_keys
            ns_nodes = Int64[]
            for (i, element_id) in enumerate(element_sets[key])
                if element_id in element_written
                    push!(ns_nodes, findfirst(x -> x == element_id, element_written))
                    continue
                end
                push!(ns_nodes, id)
                node_ids = elements[element_id]
                element_type = element_types[element_id]
                vertices = [nodes[node_id] for node_id in node_ids]
                volume = calculate_volume(string(element_type), vertices)
                center = sum(vertices) / size(vertices)[1]
                if dof == 2
                    push!(mesh_df, (x=center[1], y=center[2], volume=volume, block_id=block_id))
                else
                    push!(mesh_df, (x=center[1], y=center[2], z=center[3], volume=volume, block_id=block_id))
                end
                push!(element_written, element_id)
                id += 1
            end
            nsets[key] = ns_nodes
            block_id += 1
        end
        @info "Found $(maximum(mesh_df.block_id)) block(s)"
        @info "Found $(length(nsets)) node sets"
        @info "NodeSets: $element_sets_keys"

        mesh = nothing
        nodes = nothing
        elements = nothing
        element_sets = nothing

        return mesh_df, nsets

    elseif params["Discretization"]["Type"] == "Text File"
        return csv_reader(filename)
    else
        @error "Discretization type not supported"
        return nothing
    end

end

"""
    calculate_volume(element_type::String, vertices::Vector{Vector{Float64}})

Calculate the volume of a element.

# Arguments
- `element_type`: The element type of the element.
- `vertices`: The vertices of the element.
# Returns
- `volume`: The volume of the element.
"""
function calculate_volume(element_type::String, vertices::Vector{Vector{Float64}})
    if element_type == "Quad4"
        return area_of_polygon(vertices)
    elseif element_type == "Tet4"
        return tetrahedron_volume(vertices)
    elseif element_type == "Hex8"
        # Convert vertices to a 3x8 matrix for easier manipulation
        vertices_matrix = hcat(vertices...)
        tetrahedra_indices = delaunay(Matrix(vertices_matrix'))
        volumes = [tetrahedron_volume(tetrahedra_indices.points[tetrahedra_indices.simplices[i, :], :]) for i in 1:size(tetrahedra_indices.simplices)[1]]
        return sum(volumes)
    else
        @error "Element type $element_type currently not supported"
        return nothing
    end
end

"""
    tetrahedron_volume(vertices)

Calculate the volume of a tetrahedron.

# Arguments
- `vertices`: The vertices of the tetrahedron.
# Returns
- `volume`: The volume of the tetrahedron.
"""
function tetrahedron_volume(vertices)
    mat = hcat(vertices, ones(4))  # Augmenting matrix with ones in the fourth column
    volume = abs(det(mat) / 6)   # Using det function to calculate determinant
    return volume
end

"""
    area_of_polygon(vertices)

Calculate the area of a polygon.

# Arguments
- `vertices`: The vertices of the polygon.
# Returns
- `area`: The area of the polygon.
"""
function area_of_polygon(vertices)
    n = length(vertices)
    area = 0.0

    for i in 1:n
        j = mod(i, n) + 1
        area += (vertices[i][1] + vertices[j][1]) * (vertices[i][2] - vertices[j][2])
    end

    return abs(area) / 2.0
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
    step = params["Discretization"]["Surface Extrusion"]["Step"]
    number = params["Discretization"]["Surface Extrusion"]["Number"]

    # Finding min and max values for each dimension
    min_x, max_x = extrema(mesh.x)
    min_y, max_y = extrema(mesh.y)
    min_z = 0.0
    max_z = 0.0
    if "z" in names(mesh)
        min_z, max_z = extrema(mesh.z)
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
    volume = step * step

    id = 0

    node_sets = Dict("Extruded_1" => [], "Extruded_2" => [])

    for i in coord_max+step:step:coord_max+step*number, j in row_min:step:row_max+step#, k in min_z:step:max_z+step
        if direction == "X"
            push!(mesh, (x=i, y=j, volume=volume, block_id=block_id))
        elseif direction == "Y"
            push!(mesh, (x=j, y=i, volume=volume, block_id=block_id))
        end
        append!(node_sets["Extruded_1"], [Int64(id)])
        id += 1
    end

    block_id += 1

    for i in coord_min-step:-step:coord_min-step*number, j in row_min:step:row_max+step#, k in min_z:step:max_z+step
        if direction == "X"
            push!(mesh, (x=i, y=j, volume=volume, block_id=block_id))
        elseif direction == "Y"
            push!(mesh, (x=j, y=i, volume=volume, block_id=block_id))
        end
        append!(node_sets["Extruded_2"], [Int64(id)])
        id += 1
    end
    return mesh, node_sets
end

"""
    load_and_evaluate_mesh(params::Dict, path::String, ranksize::Int64, to::TimerOutput)

Load and evaluate the mesh data.

# Arguments
- `params::Dict`: The input parameters.
- `path::String`: The path to the mesh file.
- `ranksize::Int64`: The number of ranks.
- `to::TimerOutput`: The timer output
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
function load_and_evaluate_mesh(params::Dict, path::String, ranksize::Int64, to::TimerOutput)

    if params["Discretization"]["Type"] == "Abaqus"
        mesh, nsets = read_mesh(joinpath(path, Parameter_Handling.get_mesh_name(params)), params)
        nnodes = size(mesh, 1) + 1
        mesh, surface_ns = extrude_surface_mesh(mesh, params)
        if !isnothing(surface_ns)
            for (key, values) in surface_ns
                nsets[key] = Vector{Int64}(values .+ nnodes)
            end
        end
    else
        @debug "Read node sets"
        mesh = read_mesh(joinpath(path, Parameter_Handling.get_mesh_name(params)), params)
        nsets = get_node_sets(params, path)
    end
    check_for_duplicate_in_dataframe(mesh)
    check_types_in_dataframe(mesh)

    external_topology = nothing
    if !isnothing(get_external_topology_name(params))
        external_topology = read_external_topology(joinpath(path, get_external_topology_name(params)))
    end
    if !isnothing(external_topology)
        @info "External topology files was read."
    end
    dof::Int64 = set_dof(mesh)
    @timeit to "neighborhoodlist" nlist = create_neighborhoodlist(mesh, params, dof)
    @debug "Finished init Neighborhoodlist"
    @timeit to "apply_bond_filters" nlist, nlist_filtered_ids, bond_norm = apply_bond_filters(nlist, mesh, params, dof)
    topology = nothing
    if !isnothing(external_topology)
        @info "Create a consistent neighborhood list with external topology definition."
        nlist, topology = create_consistent_neighborhoodlist(external_topology, params["Discretization"]["Input External Topology"], nlist, dof)
    end
    @debug "Start distribution"
    if haskey(params["Discretization"], "Distribution Type")
        @timeit to "node_distribution" distribution, ptc, ntype = node_distribution(nlist, ranksize, params["Discretization"]["Distribution Type"])
    else
        @timeit to "node_distribution" distribution, ptc, ntype = node_distribution(nlist, ranksize)
    end

    el_distribution = nothing
    if haskey(params, "FEM") && !isnothing(external_topology)
        @debug "Start element distribution"
        el_distribution = element_distribution(topology, ptc, ranksize)
    end
    @debug "Finished distribution"
    @debug "Create Overlap"
    @timeit to "overlap_map" overlap_map = create_overlap_map(distribution, ptc, ranksize)
    @debug "Finished Overlap"
    return distribution, mesh, ntype, overlap_map, nlist, nlist_filtered_ids, bond_norm, dof, nsets, topology, el_distribution
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
    return neighbors(mesh, params, coor[1:dof])
end

"""
    get_number_of_neighbornodes(nlist::Vector{Vector{Int64}})

Get the number of neighbors for each node.

# Arguments
- `nlist::Vector{Vector{Int64}}`: The neighborhood list of the mesh elements.
# Returns
- `lenNlist::Vector{Int64}`: The number of neighbors for each node.
"""
function get_number_of_neighbornodes(nlist::Vector{Vector{Int64}}, filtered::Bool)
    len = length(nlist)
    lenNlist = zeros(Int64, len)
    for id in 1:len
        if !filtered && length(nlist[id]) == 0
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
                    distribution[etc[el_id]] = filter(x -> x != el_id, distribution[etc[el_id]])
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
function node_distribution(nlist::Vector{Vector{Int64}}, size::Int64, distribution_type::String="Neighbor based")

    nnodes = length(nlist)
    ntype = Dict("controllers" => Int64[], "responder" => Int64[])
    if size == 1
        distribution = [collect(1:nnodes)]
        overlap_map = [[[]]]
        ptc = Int64[]
        append!(ntype["controllers"], nnodes)
        append!(ntype["responder"], 0)
    else
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
function create_distribution_node_based(nnodes::Int64, nlist::Vector{Vector{Int64}}, size::Int64)
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
    for i in 1:size-1
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
function create_distribution_neighbor_based(nnodes::Int64, nlist::Vector{Vector{Int64}}, size::Int64)
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
function neighbors(mesh::DataFrame, params::Dict, coor::Union{Vector{Int64},Vector{String}})
    @debug "Init Neighborhoodlist"
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
    bond_intersects_disc(p0::Vector{Float64}, p1::Vector{Float64}, center::Vector{Float64}, normal::Vector{Float64}, radius::Float64)

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
function bond_intersects_disc(p0::Vector{Float64}, p1::Vector{Float64}, center::Vector{Float64}, normal::Vector{Float64}, radius::Float64)
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
    bond_intersect_infinite_plane(p0::Vector{Float64}, p1::Vector{Float64}, lower_left_corner::Vector{Float64}, normal::Vector{Float64})

Check if a line segment intersects an infinite plane.

# Arguments
- `p0::Vector{Float64}`: The start point of the line segment.
- `p1::Vector{Float64}`: The end point of the line segment.
- `lower_left_corner::Vector{Float64}`: The lower left corner of the plane.
- `normal::Vector{Float64}`: The normal of the plane.
# Returns
- `Bool`: True if the line segment intersects the plane, False otherwise.
"""
function bond_intersect_infinite_plane(p0::Vector{Float64}, p1::Vector{Float64}, lower_left_corner::Vector{Float64}, normal::Vector{Float64})
    denominator = dot((p1 - p0), normal)
    if abs(denominator) < TOLERANCE
        # Line is parallel to the plane
        # It may or may not lie on the plane
        # If it does lie on the plane, then the numerator will be zero
        # In either case, this function will return "no intersection"
        return false, undef
    end
    # The line intersects the plane

    t = dot((lower_left_corner - p0), normal) / denominator

    # Determine if the line segment intersects the plane
    if 0.0 <= t <= 1.0
        return true, p0 + t .* (p1 - p0)
    end
    # Intersection point
    return false, undef
end

"""
    bond_intersect_rectangle_plane(x::Vector{Float64}, lower_left_corner::Vector{Float64}, bottom_unit_vector::Vector{Float64}, normal::Vector{Float64}, side_length::Float64, bottom_length::Float64)

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
function bond_intersect_rectangle_plane(x::Vector{Float64}, lower_left_corner::Vector{Float64}, bottom_unit_vector::Vector{Float64}, normal::Vector{Float64}, side_length::Float64, bottom_length::Float64)
    dr::Vector{Float64} = x - lower_left_corner
    bb::Float64 = dot(dr, bottom_unit_vector)
    if 0.0 <= bb && bb / bottom_length <= 1.0
        if length(normal) == 2
            return true
        end
        ua = cross(bottom_unit_vector, normal)
        aa = dot(dr, ua)
        if 0.0 <= aa && aa / side_length <= 1.0
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
- `nlist_filtered_ids::Vector{Vector{Int64}}`: The filtered neighborhood list.
"""
function apply_bond_filters(nlist::Vector{Vector{Int64}}, mesh::DataFrame, params::Dict, dof::Int64)
    bond_filters = get_bond_filters(params)
    nlist_filtered_ids = nothing
    bond_norm = nothing
    if bond_filters[1]
        @debug "Apply bond filters"
        coor = names(mesh)[1:dof]
        nnodes = length(mesh[!, coor[1]])
        data = zeros(dof, nnodes)
        for i in 1:dof
            data[i, :] = values(mesh[!, coor[i]])
        end
        #TODO to the bottom, because right now all filters have contact if true
        contact_enabled = false
        for (name, filter) in bond_filters[2]
            contact_enabled = get(filter, "Allow Contact", false)
            if contact_enabled
                break
            end
        end
        if contact_enabled
            nlist_filtered_ids = fill(Vector{Int64}([]), nnodes)
            bond_norm = Matrix{Float64}[]
            for iID in 1:nnodes
                # bond_norm[iID] = fill(fill(1, dof), length(nlist[iID])) 
                append!(bond_norm, [Matrix{Float64}(undef, length(nlist[iID]), dof)])
                fill!(bond_norm[end], Float64(1))
            end
        end

        for (name, filter) in bond_filters[2]
            if filter["Type"] == "Disk"
                filter_flag, normal = disk_filter(nnodes, data, filter, nlist, dof)
            elseif filter["Type"] == "Rectangular_Plane"
                filter_flag, normal = rectangular_plane_filter(nnodes, data, filter, nlist, dof)
            end
            for iID in 1:nnodes
                if contact_enabled && any(x -> x == false, filter_flag[iID])
                    indices = findall(x -> x in setdiff(nlist[iID], nlist[iID][filter_flag[iID]]), nlist[iID])
                    nlist_filtered_ids[iID] = indices
                    for jID in indices
                        bond_norm[iID][jID, :] .= normal
                    end
                else
                    nlist[iID] = nlist[iID][filter_flag[iID]]
                end
            end
        end
        @debug "Finished applying bond filters"
    end
    return nlist, nlist_filtered_ids, bond_norm
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
- `filter_flag::Vector{Bool}`: The filter flag.
- `normal::Vector{Float64}`: The normal vector of the disk.
"""
function disk_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::Vector{Vector{Int64}}, dof::Int64)


    if dof == 2
        center = [filter["Center X"], filter["Center Y"]]
        normal = [filter["Normal X"], filter["Normal Y"]]
    else
        center = [filter["Center X"], filter["Center Y"], filter["Center Z"]]
        normal = [filter["Normal X"], filter["Normal Y"], filter["Normal Z"]]
    end
    #normalize vector
    normal = normal ./ norm(normal)
    filter_flag = fill(Vector{Bool}[], nnodes)
    for iID in 1:nnodes
        filter_flag[iID] = fill(true, length(nlist[iID]))
        for (jId, neighbor) in enumerate(nlist[iID])
            filter_flag[iID][jId] = !bond_intersects_disc(data[:, iID], data[:, neighbor], center, normal, filter["Radius"])
        end
    end
    return filter_flag, normal
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
- `filter_flag::Vector{Bool}`: The filter flag.
- `normal::Vector{Float64}`: The normal vector of the disk.
"""
function rectangular_plane_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::Vector{Vector{Int64}}, dof::Int64)

    if dof == 2
        normal = [filter["Normal X"], filter["Normal Y"]]
        lower_left_corner = [filter["Lower Left Corner X"], filter["Lower Left Corner Y"]]
        bottom_unit_vector = [filter["Bottom Unit Vector X"], filter["Bottom Unit Vector Y"]]
    else
        normal = [filter["Normal X"], filter["Normal Y"], filter["Normal Z"]]
        lower_left_corner = [filter["Lower Left Corner X"], filter["Lower Left Corner Y"], filter["Lower Left Corner Z"]]
        bottom_unit_vector = [filter["Bottom Unit Vector X"], filter["Bottom Unit Vector Y"], filter["Bottom Unit Vector Z"]]
    end
    #normalize vector
    normal = normal ./ norm(normal)
    bottom_unit_vector = bottom_unit_vector ./ norm(bottom_unit_vector)
    bottom_length = filter["Bottom Length"]
    side_length = filter["Side Length"]
    filter_flag::Vector{Vector{Bool}} = fill([], nnodes)
    for iID in 1:nnodes
        filter_flag[iID] = fill(true, length(nlist[iID]))
        for (jID, neighborID) in enumerate(nlist[iID])
            intersect_inf_plane, x = bond_intersect_infinite_plane(data[:, iID], data[:, neighborID], lower_left_corner, normal)
            bond_intersect = false
            if intersect_inf_plane
                bond_intersect = bond_intersect_rectangle_plane(x, lower_left_corner, bottom_unit_vector, normal, side_length, bottom_length)
            end
            filter_flag[iID][jID] = !(intersect_inf_plane && bond_intersect)
        end
    end
    return filter_flag, normal
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
