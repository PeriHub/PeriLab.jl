# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Data_manager

export check_properties
export create_bond_field
export create_constant_bond_field
export create_constant_node_field
export create_node_field
export get_all_field_keys
export get_block_list
export get_comm
export get_field
export get_local_nodes
export get_nlist
export get_nnsets
export get_nsets
export get_nnodes
export get_physics_options
export get_properties
export get_property
export get_rank
export get_nslaves
export get_max_rank
export init_property
export set_block_list
export set_glob_to_loc
export set_nmasters
export set_nset
export set_nslaves
export set_physics_options
export set_property
export set_rank
export set_max_rank
export switch_NP1_to_N
##########################
# Variables
##########################
nnodes = 0
nmasters = 0
nslaves = 0
nnsets = 0
dof = 1
block_list = Int64[]
properties = Dict{Int64,Dict{String,Any}}()
glob_to_loc = Dict{Int64,Int64}()
fields = Dict(Int64 => Dict{String,Any}(), Float64 => Dict{String,Any}(), Bool => Dict{String,Any}())
field_array_type = Dict{String,Dict{String,Any}}()
field_names = Dict{String,DataType}()
fields_to_synch = Dict{String,Any}()
nsets = Dict{String,Vector{Int}}()
overlap_map = Dict{Int64,Any}()
physicsOptions = Dict{String,Bool}("Deformed Bond Geometry" => true,
    "Deformation Gradient" => false,
    "Shape Tensor" => false,
    "Bond Associated Shape Tensor" => false,
    "Bond Associated Deformation Gradient" => false)
rank::Int64 = 0
commMPi = Any
function get_comm()
    return commMPi
end
function set_comm(value)
    global commMPi = value
end
##########################
# Material information
##########################
material_type = Dict{String,Bool}("Bond-Based" => false, "Ordinary" => false, "Correspondence" => true, "Bond-Associated" => false)
##########################

function check_property(blockId::Int64, property::String)
    if haskey(properties[blockId], property)
        return length(properties[blockId][property]) > 0
    end
    return false
end

"""
   create_bond_field(name::String, type::DataType, dof::Int64)

   Creates a bond field with the given name, data type, and degree of freedom.

   Parameters:
    - `name::String`: The name of the bond field.
    - `vartype::DataType`: The data type of the bond field.
    - `dof::Int64`: The degrees of freedom per bond.
    - `VectorOrArray::String` (optional) - Vector or Materix; Default is vector

   Returns:
   - `bond_field::Field`: The created bond field for the current time step.
   - `bond_field_np1::Field`: The created bond field for the next time step.

   Example:
   ```julia
   create_bond_field("stress", Float64, 6)  # creates a stress bond field with 6 degrees of freedom
   ```
   """

function create_bond_field(name::String, type::DataType, dof::Int64)

    return create_field(name * "N", type, "Bond_Field", dof), create_field(name * "NP1", type, "Bond_Field", dof)
end

function create_bond_field(name::String, type::DataType, VectorOrArray::String, dof::Int64)

    return create_field(name * "N", type, "Bond_Field", VectorOrArray, dof), create_field(name * "NP1", type, "Bond_Field", VectorOrArray, dof)
end

"""
   create_constant_bond_field(name::String, type::DataType, dof::Int64)

   Creates a constant bond field with the given name, data type, and degree of freedom.

   Parameters:
    - `name::String`: The name of the bond field.
    - `vartype::DataType`: The data type of the bond field.
    - `dof::Int64`: The degrees of freedom per bond.
    - `VectorOrArray::String` (optional) - Vector or Materix; Default is vector

   Returns:
   - `constant_bond_field::Field`: The created constant bond field.

   Example:
   ```julia
   create_constant_bond_field("density", Float64, 1)  # creates a density constant bond field
   ```
   """
function create_constant_bond_field(name::String, type::DataType, dof::Int64)
    return create_field(name, type, "Bond_Field", dof)
end

function create_constant_bond_field(name::String, type::DataType, VectorOrArray::String, dof::Int64)
    return create_field(name, type, "Bond_Field", VectorOrArray, dof)
end

"""
create_constant_node_field(name::String, type::DataType, dof::Int64)

Creates a constant node field with the given name, data type, and degree of freedom.

Parameters:
- `name::String`: The name of the node field.
- `vartype::DataType`: The data type of the node field.
- `dof::Int64`: The degrees of freedom per node.
- `VectorOrArray::String` (optional) - Vector or Materix; Default is vector

Returns:
- `constant_node_field::Field`: The created constant node field.

Example:
```julia
create_constant_node_field("temperature", Float64, 1)  # creates a temperature constant node field
```
"""
function create_constant_node_field(name::String, type::DataType, dof::Int64)

    return create_field(name, type, "Node_Field", dof)
end
function create_constant_node_field(name::String, type::DataType, VectorOrArray::String, dof::Int64)

    return create_field(name, type, "Node_Field", VectorOrArray, dof)
end
"""
create_field(name::String, vartype::DataType, bondNode::String, dof::Int64)

Create a field with the given `name` for the specified `vartype`. If the field already exists, return the existing field. If the field does not exist, create a new field with the specified characteristics.

# Arguments
- `name::String`: The name of the field.
- `vartype::DataType`: The data type of the field.
- `dof::Int64`: The degrees of freedom per node.

# Returns
The field with the given `name` and specified characteristics.

"""

function create_field(name::String, vartype::DataType, bondOrNode::String, dof::Int64)
    create_field(name, vartype, bondOrNode, "Vector", dof)
end

function create_field(name::String, vartype::DataType, bondOrNode::String, VectorOrArray::String, dof::Int64)
    field_dof = dof
    if haskey(fields, vartype) == false
        fields[vartype] = Dict{String,Any}()
    end
    if name in get_all_field_keys()
        return get_field(name)
    end
    if VectorOrArray == "Matrix"
        field_dof *= dof
    end
    if bondOrNode == "Node_Field"
        if dof == 1
            fields[vartype][name] = zeros(vartype, nnodes)
        else
            fields[vartype][name] = zeros(vartype, nnodes, field_dof)
        end
    elseif bondOrNode == "Bond_Field"
        nBonds = get_field("Number of Neighbors")
        if dof == 1
            fields[vartype][name] = Vector{vartype}[]
        else
            fields[vartype][name] = Matrix{vartype}[]
        end
        for i in 1:nmasters
            if dof == 1
                append!(fields[vartype][name], [Vector{vartype}(undef, nBonds[i])])
                fill!(fields[vartype][name][end], vartype(0))
            else
                append!(fields[vartype][name], [Matrix{vartype}(undef, nBonds[i], field_dof)])
                fill!(fields[vartype][name][end], vartype(0))
            end
        end
    end
    field_names[name] = vartype
    field_array_type[name] = Dict("Type" => VectorOrArray, "Dof" => dof)

    return get_field(name)
end
"""
   create_node_field(name::String, type::DataType, dof::Int64)

   Creates a node field with the given name, data type, and degree of freedom.

   Parameters:
   - `name::String`: The name of the node field.
   - `type::DataType`: The data type of the node field.
   - `dof::Int64`: The degree of freedom of each node.
   - `VectorOrArray::String` (optional) - Vector or Materix; Default is vector
   Returns:
   - `node_field::Field`: The created node field for the current time step.
   - `node_field_np1::Field`: The created node field for the next time step.

   Example:
   ```julia
   create_node_field("displacement", Float64, 3)  # creates a displacement node field with 3 degrees of freedom
   ```
   """
function create_node_field(name::String, type::DataType, dof::Int64)

    return create_field(name * "N", type, "Node_Field", dof), create_field(name * "NP1", type, "Node_Field", dof)
end

function create_node_field(name::String, type::DataType, VectorOrArray::String, dof::Int64)
    return create_field(name * "N", type, "Node_Field", VectorOrArray, dof), create_field(name * "NP1", type, "Node_Field", VectorOrArray, dof)
end

function get_all_field_keys()
    return collect(keys(field_names))
end

function get_block_list()
    return block_list
end

function get_dof()
    """
    get_dof()

    Retrieves the degree of freedom (dof) value.

    Returns:
    - `dof` (integer): The current degree of freedom value.

    Example:
    ```julia
    get_dof()  # returns the current degree of freedom
    ```
    """
    return dof
end

function get_field(name::String, time::String)

    if time == "Constant" || time == "CONSTANT"
        return get_field(name)
    end
    return get_field(name * time)

end

function get_field(name::String)

    # view() to get SubArray references
    # https://docs.julialang.org/en/v1/base/arrays/#Views-(SubArrays-and-other-view-types)
    if name in get_all_field_keys()
        if length(size(fields[field_names[name]][name])) > 1
            if field_array_type[name]["Type"] == "Matrix"
                field_dof = field_array_type[name]["Dof"]
                return view(reshape(fields[field_names[name]][name], (:, field_dof, field_dof)), :, :, :)
            end
            return view(fields[field_names[name]][name], :, :)
        end
        if field_array_type[name]["Type"] == "Matrix"
            field_dof = field_array_type[name]["Dof"]
            return view([reshape(field, (:, field_dof, field_dof)) for field in fields[field_names[name]][name]], :,)
        end
        return view(fields[field_names[name]][name], :,)
    end
    @error "Field " * name * " does not exist. Check if it is initialized as constant."
    return []
end
"""
    get_local_nodes()

    Determines the local node numbering.

    Returns:
    - `get_local_nodes` (array): returns local nodes.

    Example:
    ```julia
    get_local_nodes()  # returns local nodes or if they do not exist at the core an empty array
    ```
    """
function get_local_nodes(global_nodes)

    return [glob_to_loc[global_node] for global_node in global_nodes if global_node in keys(glob_to_loc)]

end

function get_material_type(key)
    return material_type[key]
end

function get_nlist()
    return get_field("Neighborhoodlist")
end

function get_nnodes()
    """
    get_nnodes()

    Retrieves the number of nodes.

    Returns:
    - `nmasters` (integer): The current number of nodes.

    Example:
    ```julia
    get_nmasters()  # returns the current number of nodes 
    ```
    """

    return nmasters
end

function get_NP1_to_N_Dict()
    NP1_to_N = Dict{String,String}()
    for key in get_all_field_keys()
        if occursin("NP1", key)
            NP1_to_N[key] = key[1:end-2]
        end
    end
    return NP1_to_N
end


function get_nnsets()
    """
    get_nnsets()

    Get the number of node sets.

    Returns:
    - `nnsets::Int`: The number of node sets.

    """
    return nnsets
end

"""
get_nsets()

Get the node sets

Returns:
- `nsets::Dict{String,Vector{Int64}}`: The node sets dictionary.
"""

function get_nsets()
    return nsets
end

"""
get_nslaves()

Get the the number of slave nodes

Returns:
- `nslaves::Int64`: The number of slave nodes
"""

function get_nslaves()
    return nslaves
end

function get_overlap_map()
    return overlap_map
end

function get_synch_fields()
    return fields_to_synch
end

function get_physics_options()
    if physicsOptions["Deformation Gradient"]
        physicsOptions["Shape Tensor"] = true
        physicsOptions["Deformed Bond Geometry"] = true
    end
    if physicsOptions["Bond Associated Deformation Gradient"]
        physicsOptions["Deformation Gradient"] = true
        physicsOptions["Bond Associated Shape Tensor"] = true
        physicsOptions["Shape Tensor"] = true
        physicsOptions["Deformed Bond Geometry"] = true
    end
    return physicsOptions
end

"""
    get_properties(blockId::Int64, property::Dict)

This function retrieves the value of a specified `property` for a given `blockId` if it exists in the properties dictionary.

# Arguments
- `blockId`::Int64: The identifier of the block for which to retrieve the property.
- `property`::Dict: The dictionary containing the properties for the blocks.

# Returns
- `property_value`::Any: The value associated with the specified `property` for the given `blockId`.
- `Dict()`: An empty dictionary if the specified `property` does not exist for the given `blockId`.

# Example
```julia
block_properties = Dict(
    1 => Dict("color" => "red", "size" => 10),
    2 => Dict("color" => "blue", "height" => 20)
)

# Retrieve the 'color' property for block 1
color_value = get_properties(1, block_properties, "color")  # Returns "red"

# Try to retrieve a non-existent property for block 2
non_existent_value = get_properties(2, block_properties, "width")  # Returns an empty dictionary
"""
function get_properties(blockId::Int64, property::Dict)
    if check_property(blockId, property)
        return properties[blockId][property]
    end
    return Dict()
end
"""
    get_property(blockId::Int64, property::Dict, value_name::String)

This function retrieves a specific `value_name` associated with a specified `property` for a given `blockId` if it exists in the properties dictionary.

# Arguments
- `blockId`::Int64: The identifier of the block for which to retrieve the property.
- `property`::Dict: The dictionary containing the properties for the blocks.
- `value_name`::String: The name of the value within the specified `property`.

# Returns
- `value`::Any: The value associated with the specified `value_name` within the `property` for the given `blockId`.
- `Nothing`: If the specified `blockId`, `property`, or `value_name` does not exist in the dictionary.

# Example
```julia
block_properties = Dict(
    1 => Dict("color" => Dict("value" => "red", "category" => "primary")),
    2 => Dict("color" => Dict("value" => "blue", "category" => "primary"))
)

# Retrieve the 'value' for the 'color' property of block 1
color_value = get_property(1, block_properties, "color", "value")  # Returns "red"

# Try to retrieve a non-existent value for block 2
non_existent_value = get_property(2, block_properties, "color", "width")  # Returns Nothing
"""
function get_property(blockId::Int64, property::Dict, value_name::String)
    if check_property(blockId, property)
        if value_name in keys(properties[blockId][property])
            return properties[blockId][property][value_name]
        end
    end

    return Nothing
end

"""
    get_rank()

This function returns the rank of the core.

# Returns
- `rank`::Any: The value of the `rank` variable.

# Example
```julia
current_rank = get_rank()
"""
function get_rank()
    return rank
end
"""
    get_max_rank()

This function returns the maximal rank of MPI the `max_rank`.

# Returns
- `max_rank`::Number: The value of the `max_rank` variable.

# Example
```julia
rank = get_max_rank()
"""
function get_max_rank()
    return max_rank
end

function init_property()
    for iblock in get_block_list()
        properties[iblock] = Dict{String,Dict}("Thermal Model" => Dict{String,Any}(), "Damage Model" => Dict{String,Any}(), "Material Model" => Dict{String,Any}(), "Additive Model" => Dict{String,Any}())
    end
    return collect(keys(properties[1]))
end

function set_block_list(blocks)
    global block_list = sort(unique(blocks))
end
"""
set_dof(n)

Sets the degree of freedom (dof) value globally.

Parameters:
- `n` (integer): The value to set as the degree of freedom.

Example:
```julia
set_dof(3)  # sets the degree of freedom to 3
```
"""
function set_dof(n)

    global dof = n
end
"""
set_glob_to_loc(dict)

Sets the global-to-local mapping dict globally.

Parameters:
- `dict` (array): The dict representing the global-to-local mapping.

Example:
```julia
set_glob_to_loc([1, 3, 5])  # sets the global-to-local mapping dict
```
"""
function set_glob_to_loc(dict)

    global glob_to_loc = dict
end

function set_material_type(key, value)
    if key in keys(material_type)
        global material_type[key] = value
    else
        @warn "Type " * key * " is not defined"
    end
end

"""
 set_nnodes()

 Sets the number all nodes of one core globally.

 Parameters:

 Example:
 ```
 """
function set_nnodes()
    global nnodes = nmasters + nslaves
end

"""
 set_nmasters(n::Int64)

 Sets the number of master nodes globally. For one core the number of nodes is equal to the number of master nodes.

 Parameters:
 - `n::Int64`: The value to set as the number of nodes.

 Example:
 ```julia
 set_nmasters(10)  # sets the number of nodes to 10
 ```
 """
function set_nmasters(n::Int64)
    global nmasters = n
    set_nnodes()
end
"""
set_nnsets(n::Int64)

Set the number of node sets.

Parameters:
- `n::Int64`: The number of node sets to be set.

"""
function set_nnsets(n::Int64)

    global nnsets = n
end
"""
set_nset(name, nodes)
Set the nodes associated with a named node set.

Parameters:
- `name::String`: The name of the node set.
- `nodes::Vector{Int}`: The node indices associated with the node set.

"""
function set_nset(name::String, nodes::Vector{Int})

    if name in keys(nsets)
        @warn "Node set " * name * " already defined and it is overwritten"
    end
    nsets[name] = nodes
    # set the number of node sets
    set_nnsets(length(nsets))
end

"""
set_nslaves(n::Int64)

Sets the number of slave nodes globally. For one core the number of slave is zero. Slaves hold the information of the neighbors, of one node, but are not evaluated.

Parameters:
- `n::Int64`: The value to set as the number of nodes.

Example:
```julia
set_nslaves(10)  # sets the number of slave nodes to 10
```
"""
function set_nslaves(n)
    global nslaves = n
    set_nnodes()
end

function set_overlap_map(topo)
    global overlap_map = topo
end

function set_physics_options(values)
    physicsOptions = values
end

function set_property(blockId, property, value_name, value)
    properties[blockId][property][value_name] = value
end

function set_properties(blockId, property, values)
    properties[blockId][property] = values
end

function set_rank(value)
    global rank = value
end

function set_max_rank(value)
    global max_rank = value
end

function set_synch(name, download_from_cores, upload_to_cores)

    if name in get_all_field_keys()
        field = get_field(name)
        fields_to_synch[name] = Dict{String,Any}("upload_to_cores" => upload_to_cores, "download_from_cores" => download_from_cores, "dof" => length(field[1, :]))
    elseif name * "NP1" in get_all_field_keys()
        field = get_field(name * "NP1")
        fields_to_synch[name*"NP1"] = Dict{String,Any}("upload_to_cores" => upload_to_cores, "download_from_cores" => download_from_cores, "dof" => length(field[1, :]))
    end

end

function set_fields_equal(name::String, NP1::String)
    set_fields_equal(name * NP1)
end

function set_fields_equal(NP1::String)
    NP1_to_N = get_NP1_to_N_Dict()
    if field_array_type[NP1]["Type"] == "Matrix"
        field_array_type[NP1]["Type"] = "Vector"
        field_NP1 = get_field(NP1)
        field_array_type[NP1]["Type"] = "Matrix"
        N = NP1_to_N[NP1]
        field_array_type[N]["Type"] = "Vector"
        field_N = get_field(N)
        field_array_type[N]["Type"] = "Matrix"
    else
        field_NP1 = get_field(NP1)
        N = NP1_to_N[NP1]
        field_N = get_field(N)
    end
    field_N[:] = field_NP1[:]
end

function switch_NP1_to_N()
    NP1_to_N = get_NP1_to_N_Dict()
    for NP1 in keys(NP1_to_N)
        if field_array_type[NP1]["Type"] == "Matrix"
            field_array_type[NP1]["Type"] = "Vector"
            field_NP1 = get_field(NP1)
            field_array_type[NP1]["Type"] = "Matrix"
            N = NP1_to_N[NP1]
            field_array_type[N]["Type"] = "Vector"
            field_N = get_field(N)
            field_array_type[N]["Type"] = "Matrix"
        else
            field_NP1 = get_field(NP1)
            N = NP1_to_N[NP1]
            field_N = get_field(N)
        end
        field_N[:] = field_NP1[:]
        if size(field_NP1[1]) == ()
            field_NP1[:] = fill(field_names[NP1](0), size(field_NP1))
        else
            value = 0
            if "Bond DamageNP1" == NP1
                value = 1
            end
            for fieldID in eachindex(field_NP1)
                field_NP1[fieldID] = fill(field_names[NP1](value), size(field_NP1[fieldID]))
            end
        end
    end
end
function synch_manager()
    # Liste mit den Daten die synchronisiert werden sollen -> 
    # upload; für init
    # download sum; für time int
    # down - up für bond information
end
end