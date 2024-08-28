# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Data_manager
using MPI
using DataStructures: OrderedDict

export add_active_model
export create_bond_field
export create_constant_free_size_field
export create_constant_bond_field
export create_constant_node_field
export create_node_field
export fem_active
export initialize_data
export get_active_models
export get_all_field_keys
export has_key
export get_accuracy_order
export get_block_list
export get_crit_values_matrix
export get_aniso_crit_values
export get_comm
export get_directory
export get_field
export get_field_type
export get_inverse_nlist
export get_local_nodes
export get_nlist
export get_nnsets
export get_nsets
export get_nnodes
export get_num_elements
export get_models_options
export get_properties
export get_property
export get_rank
export get_num_responder
export get_max_rank
export get_cancel
export get_damage_models
export get_material_model
export get_output_frequency
export get_rotation
export get_element_rotation
export init_property
export set_accuracy_order
export set_block_list
export set_crit_values_matrix
export set_aniso_crit_values
export set_directory
export set_inverse_nlist
export set_fem
export set_glob_to_loc
export set_damage_models
export set_material_model
export set_model_module
export set_num_controller
export set_nset
export set_num_elements
export set_num_responder
export set_models_options
export set_property
export set_rank
export set_max_rank
export set_cancel
export set_output_frequency
export set_rotation
export set_element_rotation
export switch_NP1_to_N
export synch_manager
##########################
# Variables
##########################
global active_models::OrderedDict{String,Module}
global nnodes::Int64
global num_controller::Int64
global num_responder::Int64
global num_elements::Int64
global nnsets::Int64
global dof::Int64
global fem_option::Bool
global block_list::Vector{Int64}
global distribution::Vector{Int64}
global crit_values_matrix::Array{Float64,3}
global aniso_crit_values::Dict{Int64,Vector{Float64}}
global properties::Dict{Int64,Dict{String,Any}}
global glob_to_loc::Dict{Int64,Int64}
global fields::Dict{DataType,Dict{String,Any}}
global field_array_type::Dict{String,Dict{String,Any}}
global field_types::Dict{String,DataType}
global fields_to_synch::Dict{String,Any}
global filedirectory::String
global inverse_nlist::Vector{Dict{Int64,Int64}}
global model_modules::Dict{String,Module}
global nsets::Dict{String,Vector{Int}}
global overlap_map::Dict{Int64,Any}
global models_options::Dict{String,Bool}
global output_frequency::Vector{Dict}
global accuracy_order::Int64
global rank::Int64
global commMPi::Any
global cancel::Bool
global max_rank::Int64
global silent::Bool
global rotation::Bool
global element_rotation::Bool
global damage_models::Vector{String}
global material_models::Vector{String}
##########################

"""
    initialize_data()

Initialize all parameter in the datamanager and sets them to the default values.
"""
function initialize_data()
    global nnodes
    nnodes = 0
    global num_controller
    num_controller = 0
    global num_responder
    num_responder = 0
    global num_elements
    num_elements = 0
    global nnsets
    nnsets = 0
    global dof
    dof = 2
    global fem_option
    fem_option = false
    global block_list
    block_list = []
    global distribution
    distribution = []
    global crit_values_matrix
    crit_values_matrix = fill(-1, (1, 1, 1))
    global aniso_crit_values
    aniso_crit_values = Dict()
    global properties
    properties = Dict()
    global glob_to_loc
    glob_to_loc = Dict()
    global fields
    fields = Dict(Int64 => Dict(), Float64 => Dict(), Bool => Dict())
    global field_array_type
    field_array_type = Dict()
    global field_types
    field_types = Dict()
    global fields_to_synch
    fields_to_synch = Dict()
    global filedirectory
    filedirectory = ""
    global inverse_nlist
    inverser_nlist = []
    global model_modules
    model_modules = Dict()
    global nsets
    nsets = Dict()
    global overlap_map
    overlap_map = Dict()
    global models_options
    models_options = Dict(
        "Deformed Bond Geometry" => true,
        "Deformation Gradient" => false,
        "Shape Tensor" => false,
        "Bond Associated Deformation Gradient" => false,
    )
    global output_frequency
    output_frequency = []
    global accuracy_order
    accuracy_order = 1
    global rank
    rank = 0
    global cancel
    cancel = false
    global max_rank
    max_rank = 0
    global silent
    silent = false
    global rotation
    rotation = false
    global element_rotation
    element_rotation = false
    global active_models
    active_models = Dict{String,Module}()
    global material_models = []
    global damage_models = []
end
###################################


"""
    add_active_model(module_name::Module)

Add the main modules to an array which are active.

# Arguments
- `active_module::Module`: Module of the active models.
"""
function add_active_model(key::String, active_module::Module)
    global active_models
    if !(key in keys(active_models))
        active_models[key] = active_module
    end
end


"""
    get_accuracy_order()

Returns the accuracy order for the "bond associated correspondence" implementation.

# Arguments
- `value::Int64`: The value of the accuracy_order.
"""
function get_accuracy_order()
    global accuracy_order
    return accuracy_order
end

"""
    get_comm()

Get the MPI communicator
"""
function get_comm()
    global commMPi
    return commMPi
end

function get_directory()
    global filedirectory
    return filedirectory
end

"""
    set_comm(comm::MPI.Comm)

Set the MPI communicator

# Arguments
- `comm::MPI.Comm`: MPI communicator
"""
function set_comm(comm::MPI.Comm)
    global commMPi = comm
end
##########################
# Material information
##########################
global material_type::Dict{String,Bool} = Dict(
    "Bond-Based" => false,
    "Ordinary" => false,
    "Correspondence" => true,
    "Bond-Associated" => false,
)
##########################

"""
    check_property(block_id::Int64, property::String)

Checks if the specified `property` exists for the given `block_id`.

# Arguments
- `block_id::Int64`: The ID of the block.
- `property::String`: The name of the property to check.
# Returns
- `Bool`: `true` if the property exists, `false` otherwise.
"""
function check_property(block_id::Int64, property::String)
    global properties

    if haskey(properties[block_id], property)
        return length(properties[block_id][property]) > 0
    end
    return false
end

"""
    create_bond_field(name::String, type::Type, dof::Int64)

Creates a bond field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the bond field.
- `vartype::Type`: The data type of the bond field.
- `dof::Int64`: The degrees of freedom per bond.
- `VectorOrArray::String` (optional) - Vector or Materix; Default is vector

# Returns
- `bond_field::Field`: The created bond field for the current time step.
- `bond_field_np1::Field`: The created bond field for the next time step.

Example:
```julia
create_bond_field("stress", Float64, 6)  # creates a stress bond field with 6 degrees of freedom
```
"""
function create_bond_field(
    name::String,
    type::Type,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0,
)

    return create_field(name * "N", type, "Bond_Field", dof, default_value),
    create_field(name * "NP1", type, "Bond_Field", dof, default_value)
end

function create_bond_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0,
)

    return create_field(name * "N", type, "Bond_Field", VectorOrArray, dof, default_value),
    create_field(name * "NP1", type, "Bond_Field", VectorOrArray, dof, default_value)
end

"""
    create_constant_bond_field(name::String, type::Type, dof::Int64, default_value::Union{Int64,Float64,Bool}=0))

Creates a constant bond field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the bond field.
- `vartype::Type`: The data type of the bond field.
- `dof::Int64`: The degrees of freedom per bond.
-  default_value::Union{Int64,Float64,Bool}=0) (optional) - filled with zero or false

# Returns
- `constant_bond_field::Field`: The created constant bond field.

Example:
```julia
create_constant_bond_field("density", Float64, 1)  # creates a density constant bond field
```
   """
function create_constant_bond_field(
    name::String,
    type::Type,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0,
)
    return create_field(name, type, "Bond_Field", dof, default_value)
end

function create_constant_bond_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0,
)
    return create_field(name, type, "Bond_Field", VectorOrArray, dof, default_value)
end

function create_constant_free_size_field(name::String, vartype::Type, dof::Tuple)
    return create_field(name, vartype, dof)
end

function create_free_size_field(name::String, vartype::Type, dof::Tuple)
    return create_field(name * "N", vartype, dof), create_field(name * "NP1", vartype, dof)
end

function create_field(name::String, vartype::Type, dof::Tuple)
    if !haskey(fields, vartype)
        fields[vartype] = Dict{String,Any}()
    end
    if name in get_all_field_keys()
        if size(get_field(name)) != dof
            @warn "Field $name exists already with different size. Predefined field is returned"
        end
        return get_field(name)
    end
    fields[vartype][name] = Array{vartype}(zeros(dof))
    field_types[name] = vartype
    code::String =
        "view(fields[$vartype][\"$name\"]" *
        join(repeat(", :", length(size(fields[vartype][name])))) *
        ")"
    get_function() = eval(Meta.parse(code))
    field_array_type[name] =
        Dict("Type" => "Field", "Dof" => dof, "get_function" => get_function())
    return get_field(name)
end

"""
    create_constant_node_field(name::String, type::Type, dof::Int64)

Creates a constant node field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the node field.
- `vartype::Type`: The data type of the node field.
- `dof::Int64`: The degrees of freedom per node.
- `VectorOrArray::String` (optional) - Vector or Materix; Default is vector

# Returns
- `constant_node_field::Field`: The created constant node field.

Example:
```julia
create_constant_node_field("temperature", Float64, 1)  # creates a temperature constant node field
```
"""
function create_constant_node_field(
    name::String,
    type::Type,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0,
)
    return create_field(name, type, "Node_Field", dof, default_value)
end
function create_constant_node_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0,
)
    return create_field(name, type, "Node_Field", VectorOrArray, dof, default_value)
end

"""
    create_field(name::String, vartype::Type, bondNode::String, dof::Int64, default_value::Any=0)

Create a field with the given `name` for the specified `vartype`. If the field already exists, return the existing field. If the field does not exist, create a new field with the specified characteristics.

# Arguments
- `name::String`: The name of the field.
- `vartype::Type`: The data type of the field.
- `dof::Int64`: The degrees of freedom per node.
- `default_value::Any`: The default value of the field.

# Returns
The field with the given `name` and specified characteristics.
"""
function create_field(
    name::String,
    vartype::Type,
    bondOrNode::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool},
)
    create_field(name, vartype, bondOrNode, "Vector", dof, default_value)
end

function create_field(
    name::String,
    vartype::Type,
    bondOrNode::String,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool},
)
    global nnodes
    global num_controller
    global fields
    global field_array_type
    global field_types

    field_dof = dof
    get_function = nothing

    if !haskey(fields, vartype)
        fields[vartype] = Dict{String,Any}()
    end
    if name in get_all_field_keys()
        if size(get_field(name))[1] != nnodes
            @warn "Field $name exists already with different size. Predefined field is returned"
        end
        return get_field(name)
    end
    if VectorOrArray == "Matrix"
        field_dof *= dof
    end
    if bondOrNode == "Node_Field"
        if dof == 1
            fields[vartype][name] = fill(vartype(default_value), nnodes)
            get_function = () -> view(fields[vartype][name], :)
        else
            fields[vartype][name] = fill(vartype(default_value), nnodes, field_dof)
            if VectorOrArray == "Matrix"
                get_function =
                    () -> view(reshape(fields[vartype][name], (:, dof, dof)), :, :, :)
            else
                get_function = () -> view(fields[vartype][name], :, :)
            end
        end
    elseif bondOrNode == "Bond_Field"
        nBonds = get_field("Number of Neighbors")
        if dof == 1
            fields[vartype][name] = Vector{vartype}[]
        else
            fields[vartype][name] = Matrix{vartype}[]
        end
        for i = 1:num_controller+num_responder
            if dof == 1
                append!(fields[vartype][name], [Vector{vartype}(undef, nBonds[i])])
                fill!(fields[vartype][name][end], vartype(default_value))
                get_function = () -> view(fields[vartype][name], :)
            else
                append!(
                    fields[vartype][name],
                    [Matrix{vartype}(undef, nBonds[i], field_dof)],
                )
                fill!(fields[vartype][name][end], vartype(default_value))
                if VectorOrArray == "Matrix"
                    get_function =
                        () -> view(
                            [
                                reshape(field, (:, dof, dof)) for
                                field in fields[vartype][name]
                            ],
                            :,
                        )
                else
                    get_function = () -> view(fields[vartype][name], :, :)
                end
            end
        end
    end
    field_types[name] = vartype
    field_array_type[name] =
        Dict("Type" => VectorOrArray, "Dof" => dof, "get_function" => get_function())
    return get_field(name)
end
"""
    create_node_field(name::String, type::Type, dof::Int64)

Creates a node field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the node field.
- `type::Type`: The data type of the node field.
- `dof::Int64`: The degree of freedom of each node.
- `VectorOrArray::String` (optional) - Vector or Materix; Default is vector
# Returns
- `node_field::Field`: The created node field for the current time step.
- `node_field_np1::Field`: The created node field for the next time step.

Example:
```julia
create_node_field("displacement", Float64, 3)  # creates a displacement node field with 3 degrees of freedom
```
"""
function create_node_field(name::String, type::Type, dof::Int64, default_value::Any = 0)

    return create_field(name * "N", type, "Node_Field", dof, default_value),
    create_field(name * "NP1", type, "Node_Field", dof, default_value)
end

function create_node_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Any = 0,
)
    return create_field(name * "N", type, "Node_Field", VectorOrArray, dof, default_value),
    create_field(name * "NP1", type, "Node_Field", VectorOrArray, dof, default_value)
end
"""
   fem_active()

Returns if FEM is active (true) or not (false).
"""
function fem_active()
    return fem_option
end




"""
    get_active_models()

Returns a list active model modules.
"""
function get_active_models()
    global active_models
    return active_models
end

"""
    get_all_field_keys()

Returns a list of all field keys.
"""
function get_all_field_keys()
    global field_types
    return collect(keys(field_types))
end

"""
    has_key(field_name::String)

Control if a key exists.
"""
function has_key(field_name::String)
    global field_types
    return haskey(field_types, field_name)
end

"""
    get_block_list()

Returns a list of all block IDs.
"""
function get_block_list()
    global block_list
    return block_list
end

"""
    get_crit_values_matrix()

Retrieves the critical values matrix.
"""
function get_crit_values_matrix()
    global crit_values_matrix
    return crit_values_matrix
end

"""
    get_aniso_crit_values()

Retrieves the critical values matrix.
"""
function get_aniso_crit_values()
    global aniso_crit_values
    return aniso_crit_values
end

"""
    get_dof()

Retrieves the degree of freedom (dof) value.

# Returns
- `dof` (integer): The current degree of freedom value.

Example:
```julia
get_dof()  # returns the current degree of freedom
```
"""
function get_dof()
    global dof
    return dof
end

"""
    get_field(name::String, time::String)

Returns the field with the given name and time.

# Arguments
- `name::String`: The name of the field.
- `time::String`: The time of the field.
# Returns
- `field::Field`: The field with the given name and time.
"""
function get_field(name::String, time::String, throw_error::Bool = true)

    if time == "Constant" || time == "CONSTANT"
        return get_field(name)
    end
    return get_field(name * time, throw_error)

end

"""
    get_field(name::String, throw_error::Bool=true)

Returns the field with the given name.

# Arguments
- `name::String`: The name of the field.
- `throw_error::Bool=true`: Whether to throw an error if the field does not exist.
# Returns
- `field::Field`: The field with the given name.
"""
function get_field(name::String, throw_error::Bool = true)
    global field_array_type

    if name in get_all_field_keys()
        return field_array_type[name]["get_function"]
    end
    if throw_error
        @error "Field ''" *
               name *
               "'' does not exist. Check if it is initialized as constant."
    end
    return nothing
end

"""
    get_bond_damage(time::String)
Get the bond damage

# Arguments
- `time::String`: The time of the field.
# Returns
- `damage::Field`: The bond damage field.
"""
function get_bond_damage(time::String)
    bond_damage = get_field("Bond Damage", time)
    bond_damage_aniso = get_field("Bond Damage Anisotropic", time, false)
    # return isnothing(bond_damage_aniso) ? bond_damage : bond_damage_aniso
    return bond_damage
end

"""
    get_damage(time::String)
Get the damage

# Arguments
- `time::String`: The time of the field.
# Returns
- `damage::Field`: The damage field.
"""
function get_damage(time::String)
    damage = get_field("Damage", time)
    damage_aniso = get_field("Damage Anisotropic", time, false)
    return isnothing(damage_aniso) ? damage : damage_aniso
end

"""
    get_field_type()
Get the type of a field

# Returns
- `get_field_type` (string): returns the type of a field
"""
function get_field_type(name::String)
    global field_types

    if !haskey(field_types, name)
        @error "Field ''" * name * "'' does not exist."
        return nothing
    end
    return field_types[name]
end

"""
    get_inverse_nlist()

Get the inverse of the neighborhood list.
"""
function get_inverse_nlist()
    global inverse_nlist
    return inverse_nlist
end

"""
    get_local_nodes()

Determines the local node numbering.

# Returns
- `get_local_nodes` (array): returns local nodes.

Example:
```julia
get_local_nodes()  # returns local nodes or if they do not exist at the core an empty array
```
"""
function get_local_nodes(global_nodes)
    global glob_to_loc

    return [
        glob_to_loc[global_node] for
        global_node in global_nodes if global_node in keys(glob_to_loc)
    ]

end

function get_model_module(entry::Union{String,SubString})
    global model_modules
    return model_modules[entry]
end

"""
    get_nlist()

Get the neighborhood list.
"""
function get_nlist()
    return get_field("Neighborhoodlist")
end

"""
    get_filtered_nlist()

Get the neighborhood list.
"""
function get_filtered_nlist()
    return get_field("FilteredNeighborhoodlist", false)
end

"""
    get_nnodes()

Retrieves the number of nodes.

# Returns
- `num_controller::Int64` : The current number of nodes.

Example:
```julia
get_nnodes()  # returns the current number of nodes
```
"""
function get_nnodes()
    global num_controller
    return num_controller
end

"""
    get_NP1_to_N_Dict()

Get the NP1 to N dictionary
"""
function get_NP1_to_N_Dict()
    NP1_to_N = Dict{String,String}()
    for key in get_all_field_keys()
        if occursin("NP1", key)
            NP1_to_N[key] = key[1:end-2]
        end
    end
    return NP1_to_N
end

"""
    get_nnsets()

Get the number of node sets.

# Returns
- `nnsets::Int`: The number of node sets.
"""
function get_nnsets()
    global nnsets
    return nnsets
end

"""
    get_nsets()

Get the node sets

# Returns
- `nsets::Dict{String,Vector{Int64}}`: The node sets dictionary.
"""
function get_nsets()
    global nsets
    return nsets
end

"""
    get_num_elements()

Get the the number of finite elements

# Returns
- `get_num_elements::Int64`: The number of finite elements
"""
function get_num_elements()
    global num_elements
    return num_elements
end

"""
    get_num_responder()

Get the the number of responder nodes

# Returns
- `num_responder::Int64`: The number of responder nodes
"""
function get_num_responder()
    global num_responder
    return num_responder
end

"""
    get_overlap_map()

Get the overlap map
"""
function get_overlap_map()
    global overlap_map
    return overlap_map
end

"""
    get_synch_fields()

Get the fields to synchronize
"""
function get_synch_fields()
    global fields_to_synch
    return fields_to_synch
end

"""
    get_models_options()

Get the models options
"""
function get_models_options()
    global models_options
    if models_options["Deformation Gradient"]
        models_options["Shape Tensor"] = true
        models_options["Deformed Bond Geometry"] = true
    end
    if models_options["Bond Associated Deformation Gradient"]
        models_options["Deformed Bond Geometry"] = true
    end
    return models_options
end

"""
    get_properties(block_id::Int64, property::String)

This function retrieves the value of a specified `property` for a given `block_id` if it exists in the properties dictionary.

# Arguments
- `block_id`::Int64: The identifier of the block for which to retrieve the property.
- `property`::String: The dictionary entrycontaining the properties for the blocks.

# Returns
- `property_value`::Any: The value associated with the specified `property` for the given `block_id`.
- `Dict()`: An empty dictionary if the specified `property` does not exist for the given `block_id`.

# Example
```julia
block_properties = Dict(
    1 => Dict("color" => "red", "size" => 10),
    2 => Dict("color" => "blue", "height" => 20)
)

# Retrieve the 'color' property for block 1
color_value = get_properties(1, "color")  # Returns "red"

# Try to retrieve a non-existent property for block 2
non_existent_value = get_properties(2, "width")  # Returns an empty dictionary
"""
function get_properties(block_id::Int64, property::String)
    global properties

    if check_property(block_id, property)
        return properties[block_id][property]
    end
    return Dict()
end
"""
    get_property(block_id::Int64, property::String, value_name::String)

This function retrieves a specific `value_name` associated with a specified `property` for a given `block_id` if it exists in the properties dictionary.

# Arguments
- `block_id`::Int64: The identifier of the block for which to retrieve the property.
- `property`::String: The String property type (e.g. Material model) for the blocks.
- `value_name`::String: The name of the value within the specified `property`.

# Returns
- `value`::Any: The value associated with the specified `value_name` within the `property` for the given `block_id`.
- `nothing`: If the specified `block_id`, `property`, or `value_name` does not exist in the dictionary.

# Example
```julia

"""
function get_property(block_id::Int64, property::String, value_name::String)
    global properties

    if check_property(block_id, property)
        if value_name in keys(properties[block_id][property])
            return properties[block_id][property][value_name]
        end
    end

    return nothing
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
    global rank
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
    global max_rank
    return max_rank
end

"""
    get_cancel()

This function returns the `cancel` flag.

# Returns
- `cancel`::Bool: The value of the `cancel` variable.
"""
function get_cancel()
    global cancel
    return cancel
end

"""
    get_silent()

This function returns the `silent` flag.

# Returns
- `silent`::Bool: The value of the `silent` variable.
"""
function get_silent()
    global silent
    return silent
end

"""
    get_rotation()

This function returns the `rotation` flag.

# Returns
- `rotation`::Bool: The value of the `rotation` variable.
"""
function get_rotation()
    global rotation
    return rotation
end

"""
    get_element_rotation()

This function returns the `element_rotation` flag.

# Returns
- `element_rotation`::Bool: The value of the `element_rotation` variable.
"""
function get_element_rotation()
    global element_rotation
    return element_rotation
end

"""
    get_output_frequency()

This function returns the `output_frequency` variable.

# Returns
- `output_frequency`::Any: The value of the `output_frequency` variable.
"""
function get_output_frequency()
    global output_frequency
    return output_frequency
end

"""
    get_damage_models()

This function returns the `damage_models` variable.

# Returns
- `damage_models`::Any: The value of the `damage_models` variable.
"""
function get_damage_models()
    global damage_models
    return damage_models
end

"""
    get_material_models()

This function returns the `material_models` variable.

# Returns
- `material_models`::Any: The value of the `material_models` variable.
"""
function get_material_models()
    global material_models
    return material_models
end

"""
    loc_to_glob(range::UnitRange{Int64})

Converts the local index to the global index.

# Arguments
- `range::UnitRange{Int64}`: The range of the local index.

Example:
```julia
loc_to_glob(1:10)  # converts the local index to the global index
```
"""
function loc_to_glob(range::UnitRange{Int64})
    global distribution
    return distribution[range]
end

"""
    init_property()

This function initializes the properties dictionary.

# Returns
- `keys(properties[1])`: The keys of the properties dictionary.
"""
function init_property()
    global properties

    block_list = get_block_list()
    for iblock in block_list
        properties[iblock] = Dict{String,Dict}(
            "Thermal Model" => Dict{String,Any}(),
            "Damage Model" => Dict{String,Any}(),
            "Material Model" => Dict{String,Any}(),
            "Additive Model" => Dict{String,Any}(),
            "Corrosion Model" => Dict{String,Any}(),
            "Pre Calculation Model" => Dict{String,Any}(),
        )
    end
    return collect(keys(properties[block_list[1]]))
end

"""
    set_accuracy_order(value::Int64)

Sets the accuracy order for the "bond associated correspondence" implementation.

# Arguments
- `value::Int64`: The value of the accuracy_order.
"""
function set_accuracy_order(value::Int64)
    if value < 1
        @error "Accuracy order must be greater than zero."
        return nothing
    end
    global accuracy_order = value
end


"""
    set_block_list(blocks::Union{SubArray,Vector{Int64}})

Sets the block list globally.

# Arguments
- `blocks::Union{SubArray,Vector{Int64}}`: The block list.
"""
function set_block_list(blocks::Union{SubArray,Vector{Int64}})
    global block_list = sort!(unique(blocks))
end

"""
    set_crit_values_matrix(crit_values::Array{Float64,3})

Sets the critical values matrix globally.

# Arguments
- `crit_values::Array{Float64,3}`: The critical values matrix.
"""
function set_crit_values_matrix(crit_values::Array{Float64,3})
    global crit_values_matrix = crit_values
end


"""
set_aniso_crit_values(crit_values::Dict{Int64,Any})

Sets the anisotropic critical values globally.

# Arguments
- `crit_values::Dict{Int64,Any}`: The critical values.
"""
function set_aniso_crit_values(crit_values::Dict{Int64,Any})
    global aniso_crit_values = crit_values
end

function set_directory(directory::String)
    global filedirectory = directory
end

"""
    set_distribution(values::Vector{Int64})

Sets the distribution globally.

# Arguments
- `values::Vector{Int64}`: The distribution.
"""
function set_distribution(values::Vector{Int64})
    global distribution = values
end

"""
    set_dof(n::Int64)

Sets the degree of freedom (dof) value globally.

# Arguments
- `n::Int64`: The value to set as the degree of freedom.

Example:
```julia
set_dof(3)  # sets the degree of freedom to 3
```
"""
function set_dof(n::Int64)
    global dof = n
end

"""
    set_fem(value::Bool)

Activates and deactivates the FEM option in PeriLab

# Arguments
- `value::Bool`: The value to set FEM active (true) or not (false).

Example:
```julia
set_fem(true)  # sets the fem_option to true
```
"""
function set_fem(value::Bool)
    if value
        @info "FEM is enabled"
    end
    global fem_option = value
end

"""
    set_glob_to_loc(dict)

Sets the global-to-local mapping dict globally.

# Arguments
- `dict` (array): The dict representing the global-to-local mapping.

Example:
```julia
set_glob_to_loc([1, 3, 5])  # sets the global-to-local mapping dict
```
"""
function set_glob_to_loc(dict::Dict)
    global glob_to_loc = dict
end

"""
    set_inverse_nlist(inv_nlist::Vector{Dict{Int64,Int64}})

Sets the inverse nlist globally.

# Arguments
- `inv_nlist::Vector{Dict{Int64,Int64}}`: The inverse nlist.
"""
function set_inverse_nlist(inv_nlist::Vector{Dict{Int64,Int64}})
    global inverse_nlist = inv_nlist
end

"""
    set_nnodes()

Sets the number all nodes of one core globally.

# Arguments

Example:
```
"""
function set_nnodes()
    global num_controller
    global num_responder
    global nnodes = num_controller + num_responder
end

"""
    set_num_controller(n::Int64)

Sets the number of controller nodes globally. For one core the number of nodes is equal to the number of controller nodes.

# Arguments
- `n::Int64`: The value to set as the number of nodes.

Example:
```julia
set_num_controller(10)  # sets the number of nodes to 10
```
"""
function set_num_controller(n::Int64)
    global num_controller = n
    set_nnodes()
end

"""
    set_nnsets(n::Int64)

Set the number of node sets.

# Arguments
- `n::Int64`: The number of node sets to be set.
"""
function set_nnsets(n::Int64)

    global nnsets = n
end

"""
    set_nset(name, nodes)
Set the nodes associated with a named node set.

# Arguments
- `name::String`: The name of the node set.
- `nodes::Vector{Int}`: The node indices associated with the node set.
"""
function set_nset(name::String, nodes::Vector{Int64})
    global nsets

    if name in keys(nsets)
        @warn "Node set " * name * " already defined and it is overwritten"
    end
    nsets[name] = nodes
    # set the number of node sets
    set_nnsets(length(nsets))
end

"""
    set_num_elements(n::Int64)

Sets the number of finite elements globally.

# Arguments
- `n::Int64`: The value to set as the number of finite elements.

Example:
```julia
set_num_elements(10)  # sets the number of finite elements to 10
```
"""

function set_num_elements(n::Int64)
    if n < 0
        @error "Number of elements must be positive or zero."
        return nothing
    end
    global num_elements = n
end

"""
    set_num_responder(n::Int64)

Sets the number of responder nodes globally. For one core the number of responder is zero. responder hold the information of the neighbors, of one node, but are not evaluated.

# Arguments
- `n::Int64`: The value to set as the number of nodes.

Example:
```julia
set_num_responder(10)  # sets the number of responder nodes to 10
```
"""
function set_num_responder(n::Int64)
    global num_responder = n
    set_nnodes()
end

"""
    set_overlap_map(topo)

Sets the overlap map globally.

# Arguments
- `topo`: The overlap map.
"""
function set_overlap_map(topo)
    global overlap_map = topo
end

"""
    set_models_options(values::Dict{String,Bool})

Sets the models options globally.

# Arguments
- `values::Dict{String,Bool}`: The models options.
"""
function set_models_options(values::Dict{String,Bool})
    global models_options = values
end



"""
    set_property(block_id, property, value_name, value)

Sets the value of a specified `property` for a given `block_id`.

# Arguments
- `block_id`::Int64: The identifier of the block for which to set the property.
- `property`::String: The name of the property.
- `value_name`::String: The name of the value within the specified `property`.
- `value`::Any: The value to set for the specified `value_name`.
"""
function set_property(block_id::Int64, property::String, value_name::String, value)
    global properties
    properties[block_id][property][value_name] = value
end

function set_property(property::String, value_name::String, value)
    global properties
    for block_id in eachindex(properties)
        properties[block_id][property][value_name] = value
    end
end


"""
    set_properties(block_id, property, values)

Sets the values of a specified `property` for a given `block_id`.

# Arguments
- `block_id`::Int64: The identifier of the block for which to set the property.
- `property`::String: The name of the property.
- `values`::Any: The values to set for the specified `property`.
"""
function set_properties(block_id, property, values)
    global properties
    properties[block_id][property] = values
end

"""
    set_properties(property, values)

Sets the values of a specified `property` for a all `blocks`. E.g. for FEM, because it corresponds not to a block yet,

# Arguments
- `property`::String: The name of the property.
- `values`::Any: The values to set for the specified `property`.
"""
function set_properties(property, values)
    global properties
    for id in eachindex(properties)
        properties[id][property] = values
    end
end

"""
    set_rank(value::Int64)

Sets the rank globally.

# Arguments
- `value::Int64`: The value to set as the rank.
"""
function set_rank(value::Int64)
    global rank = value
end


function set_model_module(entry::Union{String,SubString}, mod::Module)
    global model_modules[entry] = mod
end

"""
    set_max_rank(value::Int64)

Sets the maximum rank globally.

# Arguments
- `value::Int64`: The value to set as the maximum rank.
"""
function set_max_rank(value::Int64)
    global max_rank = value
end

"""
    set_cancel(value::Int64)

Sets the cancel flag.

# Arguments
- `value::Bool`: The cancel flag.
"""
function set_cancel(value::Bool)
    global cancel = value
end

"""
    set_silent(value::Int64)

Sets the silent flag.

# Arguments
- `value::Bool`: The silent flag.
"""
function set_silent(value::Bool)
    global silent = value
end

"""
    set_rotation(value::Int64)

Sets the rotation flag.

# Arguments
- `value::Bool`: The rotation flag.
"""
function set_rotation(value::Bool)
    global rotation = value
end

"""
    set_element_rotation(value::Int64)

Sets the element_rotation flag.

# Arguments
- `value::Bool`: The element_rotation flag.
"""
function set_element_rotation(value::Bool)
    global element_rotation = value
end

"""
    set_output_frequency(value)

Sets the output frequency globally.

# Arguments
- `value`: The value to set as the output frequency.
"""
function set_output_frequency(value)
    global output_frequency = value
end

"""
    set_damage_models(value)

Sets the damage models globally.

# Arguments
- `value`: The value to set as the damage models.
"""
function set_damage_models(value)
    global damage_models
    if !(value in damage_models)
        push!(damage_models, value)
    end
end

"""
    set_material_models(value)

Sets the material models globally.

# Arguments
- `value`: The value to set as the material models.
"""
function set_material_models(value)
    global material_models
    if !(value in material_models)
        push!(material_models, value)
    end
end

"""
    set_synch(name, download_from_cores, upload_to_cores)

Sets the synchronization dictionary globally.

# Arguments
- `name`::String: The name of the field.
- `download_from_cores`::Bool: Whether to download the field from the cores.
- `upload_to_cores`::Bool: Whether to upload the field to the cores.
"""
function set_synch(name, download_from_cores, upload_to_cores)
    global fields_to_synch

    if name in get_all_field_keys()
        field = get_field(name)
        fields_to_synch[name] = Dict{String,Any}(
            "upload_to_cores" => upload_to_cores,
            "download_from_cores" => download_from_cores,
            "dof" => length(field[1, :]),
        )
    elseif name * "NP1" in get_all_field_keys()
        field = get_field(name * "NP1")
        fields_to_synch[name*"NP1"] = Dict{String,Any}(
            "upload_to_cores" => upload_to_cores,
            "download_from_cores" => download_from_cores,
            "dof" => length(field[1, :]),
        )
    end

end

"""
    switch_NP1_to_N()

Switches the fields from NP1 to N.
"""
function switch_NP1_to_N()
    global field_types
    global field_array_type

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

        if size(field_NP1[1]) == () # vector

            copyto!(field_N, field_NP1)
            fill!(field_NP1, field_types[NP1](0))
        else # matrix
            value = 0
            for fieldID in eachindex(field_NP1)
                copyto!(field_N[fieldID], field_NP1[fieldID])
                if "Bond DamageNP1" != NP1
                    # value = 1
                    fill!(field_NP1[fieldID], field_types[NP1](value))
                end
            end
        end
    end
end

"""
    synch_manager(synchronise_field, direction::String)

Synchronises the fields.

# Arguments
- `synchronise_field`: The function to synchronise the field.
- `direction::String`: The direction of the synchronisation.
"""
function synch_manager(synchronise_field, direction::String)
    global overlap_map

    synch_fields = get_synch_fields()
    for synch_field in keys(synch_fields)
        synchronise_field(
            get_comm(),
            synch_fields,
            overlap_map,
            get_field,
            synch_field,
            direction,
        )
    end
end
end
