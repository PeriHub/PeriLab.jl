# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Data_manager
using MPI
using DataStructures: OrderedDict
include("../Support/Helpers.jl")
using .Helpers: fill_in_place!

export add_active_model
export create_bond_field
export create_constant_free_size_field
export create_constant_bond_field
export create_constant_node_field
export create_constant_element_field
export create_node_field
export fem_active
export initialize_data
export get_active_models
export get_all_field_keys
export has_key
export get_bc_free_dof
export get_step
export get_iteration
export get_accuracy_order
export get_aniso_crit_values
export get_block_list
export get_crit_values_matrix
export get_comm
export get_coupling_dict
export get_coupling_fe_nodes
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
export get_pre_calculation_order
export get_properties
export get_property
export get_rank
export get_num_responder
export get_max_rank
export get_cancel
export get_damage_models
export get_output_frequency
export get_rotation
export get_element_rotation
export init_properties
export remove_active_model
export set_step
export set_iteration
export set_accuracy_order
export set_bc_free_dof
export set_block_list
export set_crit_values_matrix
export set_coupling_dict
export set_coupling_fe_nodes
export set_aniso_crit_values
export set_directory
export set_inverse_nlist
export set_fem
export set_glob_to_loc
export set_damage_models
export set_model_module
export set_num_controller
export set_nset
export set_num_elements
export set_num_responder
export set_pre_calculation_order
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
const fields = Dict()
const data = Dict()
##########################

"""
    initialize_data()

Initialize all parameter in the datamanager and sets them to the default values.
"""
function initialize_data()
    data["step"] = 0
    data["max_step"] = 0
    data["nnodes"] = 0
    data["num_controller"] = 0
    data["num_responder"] = 0
    data["num_elements"] = 0
    data["nnsets"] = 0
    data["dof"] = 2
    data["fem_option"] = false
    data["block_list"] = Vector{String}()
    data["distribution"] = []
    data["crit_values_matrix"] = fill(-1, (1, 1, 1))
    data["aniso_crit_values"] = Dict()
    data["properties"] = OrderedDict()
    data["glob_to_loc"] = Dict()
    data["field_array_type"] = Dict()
    data["field_types"] = Dict()
    data["field_names"] = Vector{String}([])
    data["fields_to_synch"] = Dict()
    data["local_fields_to_synch"] = Dict(
        "Material Model" => Dict(),
        "Damage Model" => Dict(),
        "Thermal Model" => Dict(),
        "Pre Calculation Model" => Dict(),
        "Corrosion Model" => Dict(),
        "Additive Model" => Dict(),
    )
    data["filedirectory"] = ""
    data["inverse_nlist"] = []
    data["model_modules"] = OrderedDict{String,Module}()
    data["nsets"] = Dict{String,Vector{Int64}}()
    data["overlap_map"] = Dict()
    data["pre_calculation_order"] = [
        "Deformed Bond Geometry",
        "Shape Tensor",
        "Deformation Gradient",
        "Bond Associated Correspondence",
    ]
    data["coupling_dict"] = Dict{Int64,Int64}()
    data["output_frequency"] = []
    data["accuracy_order"] = 1
    data["rank"] = 0
    data["cancel"] = false
    data["max_rank"] = 0
    data["silent"] = false
    data["rotation"] = false
    data["element_rotation"] = false
    data["active_models"] = OrderedDict{String,Module}()
    data["all_active_models"] = OrderedDict{String,Module}()
    data["material_models"] = []
    data["damage_models"] = []
    data["NP1_to_N"] = Dict{String,Vector{}}()
    data["coupling_fe_nodes"] = []
    data["BC_free_dof"] = []
    fields[Int64] = Dict()
    fields[Float64] = Dict()
    fields[Bool] = Dict()
end
###################################

"""
    set_step(step::Int64)

Set the step of the simulation.

# Arguments
- `step::Int64`: The step of the simulation.

"""
function set_step(step::Union{Int64,Nothing})
    data["step"] = step
end

"""
    get_step()

Get the step of the simulation.

# Returns
- `Int64`: The step of the simulation.

"""
function get_step()
    return data["step"]
end

"""
    set_max_step(max_step::Int64)

Set the max_step of the simulation.

# Arguments
- `max_step::Int64`: The max_step of the simulation.

"""
function set_max_step(max_step::Union{Int64,Nothing})
    data["max_step"] = max_step
end
"""
    get_bc_free_dof()

Get all dof without displacment boundary conditions.

# Returns
- `Vector{Tuple{Int64, Int64}}`: The point and dof without boundary condition.

"""
function get_bc_free_dof()
    return data["BC_free_dof"]
end

"""
    set_bc_free_dof(values::Vector{Tuple{Int64, Int64}})

Set all dof without displacment boundary conditions.

# Returns
-

"""
function set_bc_free_dof(values::Vector{Int64})
    data["BC_free_dof"] = values
end

"""
    get_max_step()

Get the max_step of the simulation.

# Returns
- `Int64`: The max_step of the simulation.

"""
function get_max_step()
    return data["max_step"]
end

"""
    set_iteration(iteration::Int64)

Set the iteration of the simulation.

# Arguments
- `iteration::Int64`: The iteration of the simulation.

"""
function set_iteration(iteration::Int64)
    data["iteration"] = iteration
end

"""
    get_step()

Get the iteration of the simulation.

# Returns
- `Int64`: The iteration of the simulation.

"""
function get_iteration()
    return data["iteration"]
end

"""
    add_active_model(key::String, module_name::Module)

Add the main modules to an OrderedDict which are active.

# Arguments
- `key::String`: Name of the model.
- `active_module::Module`: Module of the active models.
"""
function add_active_model(key::String, active_module::Module, all::Bool = false)
    if all
        if !(key in keys(data["all_active_models"]))
            data["all_active_models"][key] = active_module
        end
    else
        if !(key in keys(data["active_models"]))
            data["active_models"][key] = active_module
        end
    end
end


"""
    get_accuracy_order()

Returns the accuracy order for the "bond associated correspondence" implementation.

# Arguments
- `value::Int64`: The value of the accuracy_order.
"""
function get_accuracy_order()
    return data["accuracy_order"]
end

"""
    get_comm()

Get the MPI communicator
"""
function get_comm()
    return data["commMPi"]
end


"""
    get_coupling_dict()

Get the PD - FE coupling dict
"""
function get_coupling_dict()
    return data["coupling_dict"]
end

"""
    get_coupling_fe_nodes()

Get the FE nodes involved in the coupling
"""
function get_coupling_fe_nodes()
    return data["coupling_fe_nodes"]
end




function get_directory()
    return data["filedirectory"]
end

"""
    set_comm(comm::MPI.Comm)

Set the MPI communicator

# Arguments
- `comm::MPI.Comm`: MPI communicator
"""
function set_comm(comm::MPI.Comm)
    data["commMPi"] = comm
end

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
    haskey(data["properties"][block_id], property) &&
        !isempty(data["properties"][block_id][property])
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
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    data["NP1_to_N"][name] = [name * "N", name * "NP1", zero(type)]
    return create_field(name * "N", type, "Bond_Field", dof, default_value),
    create_field(name * "NP1", type, "Bond_Field", dof, default_value)
end

function create_bond_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)

    data["NP1_to_N"][name] = [name * "N", name * "NP1", zero(type)]
    return create_field(name * "N", type, "Bond_Field", dof, default_value, VectorOrArray),
    create_field(name * "NP1", type, "Bond_Field", dof, default_value, VectorOrArray)
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
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    return create_field(name, type, "Bond_Field", dof, default_value)
end

function create_constant_bond_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    return create_field(name, type, "Bond_Field", dof, default_value, VectorOrArray)
end

function create_constant_free_size_field(
    name::String,
    vartype::Type,
    dof::Tuple,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    return create_field(name, vartype, "Free_Size_Field", dof, default_value)
end

function create_free_size_field(
    name::String,
    vartype::Type,
    dof::Tuple,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    data["NP1_to_N"][name] = [name * "N", name * "NP1", zero(vartype)]
    return create_field(name * "N", vartype, "Free_Size_Field", dof, default_value),
    create_field(name * "NP1", vartype, "Free_Size_Field", dof, default_value)
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
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    return create_field(name, type, "Node_Field", dof, default_value)
end
function create_constant_node_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    return create_field(name, type, "Node_Field", dof, default_value, VectorOrArray)
end

"""
    create_constant_element_field(name::String, type::Type, dof::Int64)

Creates a constant element field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the element field.
- `vartype::Type`: The data type of the element field.
- `dof::Int64`: The degrees of freedom per element.
- `VectorOrArray::String` (optional) - Vector or Materix; Default is vector

# Returns
- `constant_element_field::Field`: The created constant element field.

Example:
```julia
create_constant_element_field("temperature", Float64, 1)  # creates a temperature constant element field
```
"""
function create_constant_element_field(
    name::String,
    type::Type,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    return create_field(name, type, "Element_Field", dof, default_value)
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
function create_node_field(
    name::String,
    type::Type,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)

    data["NP1_to_N"][name] = [name * "N", name * "NP1", zero(type)]
    return create_field(name * "N", type, "Node_Field", dof, default_value),
    create_field(name * "NP1", type, "Node_Field", dof, default_value)
end

function create_node_field(
    name::String,
    type::Type,
    VectorOrArray::String,
    dof::Int64,
    default_value::Union{Int64,Float64,Bool} = 0.0,
)
    data["NP1_to_N"][name] = [name * "N", name * "NP1", zero(type)]
    return create_field(name * "N", type, "Node_Field", dof, default_value, VectorOrArray),
    create_field(name * "NP1", type, "Node_Field", dof, default_value, VectorOrArray)
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
    bond_or_node::String,
    dof::Union{Int64,Tuple},
    default_value::Union{Int64,Float64,Bool},
    VectorOrArray::String = "Vector",
)
    if has_key(name)
        if size(_get_field(name), 1) != data["nnodes"]
            @warn "Field $name exists already with different size. Predefined field is returned"
        end
        return _get_field(name)
    end
    value = vartype(default_value)
    if bond_or_node == "Node_Field"
        if dof == 1
            fields[vartype][name] = fill(value, data["nnodes"])
        else
            if VectorOrArray == "Matrix"
                fields[vartype][name] = fill(value, (data["nnodes"], dof, dof))
            else
                fields[vartype][name] = fill(value, data["nnodes"], dof)
                # fields[vartype][name] = [fill(value,dof) for j=1:data["nnodes"]]
            end
        end
    elseif bond_or_node == "Bond_Field"
        nBonds = _get_field("Number of Neighbors")
        if dof == 1
            fields[vartype][name] = [fill(value, n) for n in nBonds]
        else
            if VectorOrArray == "Matrix"
                fields[vartype][name] = [fill(value, (n, dof, dof)) for n in nBonds]
            else
                # fields[vartype][name] = [fill(value, (n, dof)) for n in nBonds]
                fields[vartype][name] = [[fill(value, dof) for j = 1:n] for n in nBonds]
            end
        end
    elseif bond_or_node == "Element_Field"
        nElements = _get_field("Number of Element Neighbors")
        if dof == 1
            fields[vartype][name] = [fill(value, n) for n in nElements]
        else
            fields[vartype][name] = [fill(value, (n, dof)) for n in nElements]
        end
    elseif bond_or_node == "Free_Size_Field"
        fields[vartype][name] = Array{vartype}(zeros(dof))
    end
    get_function = () -> fields[vartype][name]
    data["field_types"][name] = vartype
    data["field_array_type"][name] =
        Dict("Type" => VectorOrArray, "Dof" => dof, "get_function" => get_function)
    data["field_names"] = Vector{String}(collect(keys(data["field_types"])))
    return get_function()
end

"""
fem_active()

Returns if FEM is active (true) or not (false).
"""
function fem_active()
    return data["fem_option"]
end

"""
    get_active_models()

Returns a list active model modules.
"""
function get_active_models(all::Bool = false)
    return all ? data["all_active_models"] : data["active_models"]
end

"""
    get_all_field_keys()

Returns a list of all field keys.
"""
function get_all_field_keys()
    return data["field_names"]
end

"""
    has_key(field_name::String)

Control if a key exists.
"""
function has_key(field_name::String)
    return field_name in data["field_names"]
end

"""
    get_block_list()

Returns a list of all block IDs.
"""
function get_block_list()
    return data["block_list"]
end

"""
    get_crit_values_matrix()

Retrieves the critical values matrix.
"""
function get_crit_values_matrix()
    return data["crit_values_matrix"]
end

"""
    get_aniso_crit_values()

Retrieves the critical values matrix.
"""
function get_aniso_crit_values()
    return data["aniso_crit_values"]
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
    return data["dof"]
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
function get_field(name::String, time::String = "constant")

    if time == "constant"
        return _get_field(name)
    elseif time == "N"
        try
            return _get_field(data["NP1_to_N"][name][1])
        catch
            @error "Field ''" *
                   name *
                   "'' does not exist. Check if it is initialized as constant."
            return nothing
        end
    elseif time == "NP1"
        try
            return _get_field(data["NP1_to_N"][name][2])
        catch
            @error "Field ''" *
                   name *
                   "'' does not exist. Check if it is initialized as constant."
            return nothing
        end
    else
        @error "Time $time is not supported. Use 'constant', 'N', or 'NP1'"
        return nothing
    end
end

"""
    get_field_if_exists(name::String, time::String)

Returns the field with the given name if it exists.

# Arguments
- `name::String`: The name of the field.
- `time::String`: The time of the field.
# Returns
- `field::Field`: The field with the given name and time.
"""
function get_field_if_exists(name::String, time::String = "constant")
    return has_key(name) ? get_field(name, time) : nothing
end
"""
    _get_field(name::String)

Returns the field with the given name.

# Arguments
- `name::String`: The name of the field.
# Returns
- `field::Field`: The field with the given name.
"""
function _get_field(name::String)
    try
        return data["field_array_type"][name]["get_function"]()
    catch
        @error "Field ''" *
               name *
               "'' does not exist. Check if it is initialized as non-constant."
        return nothing
    end
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
    # bond_damage_aniso = get_field("Bond Damage Anisotropic", time, false)
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
    if "Damage Anisotropic" in get_all_field_keys()
        damage_aniso = get_field("Damage Anisotropic", time)
        return damage_aniso
    end
    damage = get_field("Damage", time)
    return damage
end

"""
    get_field_type()
Get the type of a field

# Returns
- `get_field_type` (string): returns the type of a field
"""
function get_field_type(name::String)
    if name in data["field_names"]
        return data["field_types"][name]
    end
    @error "Field ''" * name * "'' does not exist."
    return nothing
end

"""
    get_inverse_nlist()

Get the inverse of the neighborhood list.
"""
function get_inverse_nlist()
    return data["inverse_nlist"]
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
    return [
        data["glob_to_loc"][global_node] for
        global_node in global_nodes if global_node in keys(data["glob_to_loc"])
    ]

end

function get_model_module(entry::Union{String,SubString})
    return data["model_modules"][entry]
end

"""
    get_nlist()

Get the neighborhood list.
"""
function get_nlist()
    return _get_field("Neighborhoodlist")
end

"""
    get_filtered_nlist()

Get the neighborhood list.
"""
function get_filtered_nlist()
    if !has_key("FilteredNeighborhoodlist")
        return nothing
    end
    return _get_field("FilteredNeighborhoodlist")
end

"""
    get_nnodes()

Retrieves the number of nodes.

# Returns
- `num_controller::Int64` : The current number of nodes.

Example:
```julia
get_nnodes()  # returns the current number of controler nodes. The neighbors are not included
```
"""
function get_nnodes()
    return data["num_controller"]
end

"""
    get_nnsets()

Get the number of node sets.

# Returns
- `nnsets::Int`: The number of node sets.
"""
function get_nnsets()
    return data["nnsets"]
end

"""
    get_nsets()

Get the node sets

# Returns
- `nsets::Dict{String,Vector{Int64}}`: The node sets dictionary.
"""
function get_nsets()
    return data["nsets"]
end

"""
    get_num_elements()

Get the the number of finite elements

# Returns
- `get_num_elements::Int64`: The number of finite elements
"""
function get_num_elements()
    return data["num_elements"]
end

"""
    get_num_responder()

Get the the number of responder nodes

# Returns
- `num_responder::Int64`: The number of responder nodes
"""
function get_num_responder()
    return data["num_responder"]
end

"""
    get_overlap_map()

Get the overlap map
"""
function get_overlap_map()
    return data["overlap_map"]
end

"""
    get_synch_fields()

Get the fields to synchronize
"""
function get_synch_fields()
    return data["fields_to_synch"]
end

"""
    get_local_synch_fields(model::String)

    model - class of models; before computation of these models the synchronisation occurs
    Get the fields to synchronize
"""
function get_local_synch_fields(model::String)
    return data["local_fields_to_synch"][model]
end

"""
    get_pre_calculation_order()

return the order of the pre calculation.
"""
function get_pre_calculation_order()
    return data["pre_calculation_order"]
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
    if check_property(block_id, property)
        return data["properties"][block_id][property]
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
    if check_property(block_id, property) &&
       haskey(data["properties"][block_id][property], value_name)
        return data["properties"][block_id][property][value_name]
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
    return data["rank"]
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
    return data["max_rank"]
end

"""
    get_cancel()

This function returns the `cancel` flag.

# Returns
- `cancel`::Bool: The value of the `cancel` variable.
"""
function get_cancel()
    return data["cancel"]
end

"""
    get_silent()

This function returns the `silent` flag.

# Returns
- `silent`::Bool: The value of the `silent` variable.
"""
function get_silent()
    return data["silent"]
end

"""
    get_rotation()

This function returns the `rotation` flag.

# Returns
- `rotation`::Bool: The value of the `rotation` variable.
"""
function get_rotation()
    return data["rotation"]
end

"""
    get_element_rotation()

This function returns the `element_rotation` flag.

# Returns
- `element_rotation`::Bool: The value of the `element_rotation` variable.
"""
function get_element_rotation()
    return data["element_rotation"]
end

"""
    get_output_frequency()

This function returns the `output_frequency` variable.

# Returns
- `output_frequency`::Any: The value of the `output_frequency` variable.
"""
function get_output_frequency()
    return data["output_frequency"]
end

"""
    get_damage_models()

This function returns the `damage_models` variable.

# Returns
- `damage_models`::Any: The value of the `damage_models` variable.
"""
function get_damage_models()
    return data["damage_models"]
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
    return data["distribution"][range]
end

"""
    init_properties()

This function initializes the properties dictionary. Order of dictionary defines, in which order the models are called later on.

# Returns
- `keys(properties[1])`: The keys of the properties dictionary in defined order for the Model_Factory.jl.
"""
function init_properties()

    block_list = get_block_list()
    for iblock = 1:length(block_list)
        data["properties"][iblock] = OrderedDict{String,Dict}(
            "Additive Model" => Dict{String,Any}(),
            "Damage Model" => Dict{String,Any}(),
            "Pre Calculation Model" => Dict{String,Any}(),
            "Thermal Model" => Dict{String,Any}(),
            "Corrosion Model" => Dict{String,Any}(),
            "Material Model" => Dict{String,Any}(),
        )
    end
    return collect(keys(data["properties"][1]))
end

"""
    remove_active_model(module_name::Module)

Removes main modules from OrderedDict.

# Arguments
- `key::String`: Key of the entry.
"""
function remove_active_model(key::String)
    delete!(data["active_models"], key)
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
    data["accuracy_order"] = value
end


"""
    set_block_list(blocks::Union{SubArray,Vector{Int64}})

Sets the block list globally.

# Arguments
- `blocks::Union{SubArray,Vector{Int64}}`: The block list.
"""
function set_block_list(blocks::Vector{String})
    data["block_list"] = blocks
end

"""
    set_crit_values_matrix(crit_values::Array{Float64,3})

Sets the critical values matrix globally.

# Arguments
- `crit_values::Array{Float64,3}`: The critical values matrix.
"""
function set_crit_values_matrix(crit_values::Array{Float64,3})
    data["crit_values_matrix"] = crit_values
end


"""
    set_coupling_dict(coupling_dict::Dict{Int64,Int64})

Sets the FE - PD couplings. PD nodes -> FE Elements.

# Arguments
- `coupling_dict::Dict{Int64,Int64}`: The coupling dictionary.
"""
function set_coupling_dict(coupling_dict::Dict{Int64,Int64})
    data["coupling_dict"] = coupling_dict
end


"""
    set_coupling_fe_nodes()

Get the FE nodes involved in the coupling
"""
function set_coupling_fe_nodes(values::Vector{Int64})
    data["coupling_fe_nodes"] = values
end

"""
set_aniso_crit_values(crit_values::Dict{Int64,Any})

Sets the anisotropic critical values globally.

# Arguments
- `crit_values::Dict{Int64,Any}`: The critical values.
"""
function set_aniso_crit_values(crit_values::Dict{Int64,Any})
    data["aniso_crit_values"] = crit_values
end

function set_directory(directory::String)
    data["filedirectory"] = directory
end

"""
    set_distribution(values::Vector{Int64})

Sets the distribution globally.

# Arguments
- `values::Vector{Int64}`: The distribution.
"""
function set_distribution(values::Vector{Int64})
    data["distribution"] = values
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
    data["dof"] = n
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
    data["fem_option"] = value
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
    data["glob_to_loc"] = dict
end

"""
    set_inverse_nlist(inv_nlist::Vector{Dict{Int64,Int64}})

Sets the inverse nlist globally.

# Arguments
- `inv_nlist::Vector{Dict{Int64,Int64}}`: The inverse nlist.
"""
function set_inverse_nlist(inv_nlist::Vector{Dict{Int64,Int64}})
    data["inverse_nlist"] = inv_nlist
end

"""
    set_nnodes()

Sets the number all nodes of one core globally.

# Arguments

Example:
```
"""
function set_nnodes()
    data["num_controller"]
    data["num_responder"]
    data["nnodes"] = data["num_controller"] + data["num_responder"]
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
    data["num_controller"] = n
    set_nnodes()
end

"""
    set_nnsets(n::Int64)

Set the number of node sets.

# Arguments
- `n::Int64`: The number of node sets to be set.
"""
function set_nnsets(n::Int64)

    data["nnsets"] = n
end

"""
    set_nset(name, nodes)
Set the nodes associated with a named node set.

# Arguments
- `name::String`: The name of the node set.
- `nodes::Vector{Int}`: The node indices associated with the node set.
"""
function set_nset(name::String, nodes::Vector{Int64})
    if name in keys(data["nsets"])
        @warn "Node set " * name * " already defined and it is overwritten"
    end
    data["nsets"][name] = nodes
    # set the number of node sets
    set_nnsets(length(data["nsets"]))
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
    data["num_elements"] = n
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
    data["num_responder"] = n
    set_nnodes()
end

"""
    set_overlap_map(topo)

Sets the overlap map globally.

# Arguments
- `topo`: The overlap map.
"""
function set_overlap_map(topo)
    data["overlap_map"] = topo
end

"""
    set_pre_calculation_order(values::Vector{String})

Sets the order of the pre calculation options globally.

# Arguments
- `values::Vector{String}`: The order of models.
"""
function set_pre_calculation_order(values::Vector{String})
    data["pre_calculation_order"] = values
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
    data["properties"][block_id][property][value_name] = value
end

function set_property(property::String, value_name::String, value)
    for block_id in eachindex(data["properties"])
        data["properties"][block_id][property][value_name] = value
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
    data["properties"][block_id][property] = values
end

"""
    set_properties(property, values)

Sets the values of a specified `property` for a all `blocks`. E.g. for FEM, because it corresponds not to a block yet,

# Arguments
- `property`::String: The name of the property.
- `values`::Any: The values to set for the specified `property`.
"""
function set_properties(property, values)
    for id in eachindex(data["properties"])
        data["properties"][id][property] = values
    end
end

"""
    set_rank(value::Int64)

Sets the rank globally.

# Arguments
- `value::Int64`: The value to set as the rank.
"""
function set_rank(value::Int64)
    data["rank"] = value
end


function set_model_module(entry::Union{String,SubString}, mod::Module)
    data["model_modules"][entry] = mod
end

"""
    set_max_rank(value::Int64)

Sets the maximum rank globally.

# Arguments
- `value::Int64`: The value to set as the maximum rank.
"""
function set_max_rank(value::Int64)
    data["max_rank"] = value
end

"""
    set_cancel(value::Int64)

Sets the cancel flag.

# Arguments
- `value::Bool`: The cancel flag.
"""
function set_cancel(value::Bool)
    data["cancel"] = value
end

"""
    set_silent(value::Int64)

Sets the silent flag.

# Arguments
- `value::Bool`: The silent flag.
"""
function set_silent(value::Bool)
    data["silent"] = value
end

"""
    set_rotation(value::Int64)

Sets the rotation flag.

# Arguments
- `value::Bool`: The rotation flag.
"""
function set_rotation(value::Bool)
    data["rotation"] = value
end

"""
    set_element_rotation(value::Int64)

Sets the element_rotation flag.

# Arguments
- `value::Bool`: The element_rotation flag.
"""
function set_element_rotation(value::Bool)
    data["element_rotation"] = value
end

"""
    set_output_frequency(value)

Sets the output frequency globally.

# Arguments
- `value`: The value to set as the output frequency.
"""
function set_output_frequency(value)
    data["output_frequency"] = value
end

"""
    set_damage_models(value)

Sets the damage models globally.

# Arguments
- `value`: The value to set as the damage models.
"""
function set_damage_models(value)
    if !(value in data["damage_models"])
        push!(data["damage_models"], value)
    end
end

"""
    set_material_models(value)

Sets the material models globally.

# Arguments
- `value`: The value to set as the material models.
"""
function set_material_models(value)
    if !(value in data["material_models"])
        push!(data["material_models"], value)
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
function set_synch(name, download_from_cores, upload_to_cores, dof = 0)
    if name in get_all_field_keys()
        field = _get_field(name)
        if dof == 0
            dof = length(field[1, :, :])
        end
        data["fields_to_synch"][name] = Dict{String,Any}(
            "upload_to_cores" => upload_to_cores,
            "download_from_cores" => download_from_cores,
            "dof" => dof,
            "time" => "constant",
        )
    elseif name * "NP1" in get_all_field_keys()
        field = get_field(name, "NP1")
        if dof == 0
            dof = length(field[1, :, :])
        end
        data["fields_to_synch"][name] = Dict{String,Any}(
            "upload_to_cores" => upload_to_cores,
            "download_from_cores" => download_from_cores,
            "dof" => length(field[1, :, :]),
            "time" => "NP1",
        )
    end

end


"""
    set_local_synch(model, name, download_from_cores, upload_to_cores, dof=0)

Sets the synchronization dictionary locally during the model update process. Should be used carefully, to avoid unessary communication.

# Arguments
- `name`::String: The name of the field.
- `download_from_cores`::Bool: Whether to download the field from the cores.
- `upload_to_cores`::Bool: Whether to upload the field to the cores.
"""
function set_local_synch(model, name, download_from_cores, upload_to_cores, dof = 0)
    if !haskey(data["local_fields_to_synch"], model)
        @error "Model $model is not defined. If it is a new model type, please add it to data[\"local_fields_to_synch\"] in the datamanager."
        return nothing
    end
    if name in get_all_field_keys()
        field = _get_field(name)
        if dof == 0
            dof = length(field[1, :, :])
        end
        data["local_fields_to_synch"][model][name] = Dict{String,Any}(
            "upload_to_cores" => upload_to_cores,
            "download_from_cores" => download_from_cores,
            "dof" => dof,
            "time" => "constant",
        )
    elseif name * "NP1" in get_all_field_keys()
        field = get_field(name, "NP1")
        if dof == 0
            dof = length(field[1, :, :])
        end
        data["local_fields_to_synch"][model][name] = Dict{String,Any}(
            "upload_to_cores" => upload_to_cores,
            "download_from_cores" => download_from_cores,
            "dof" => length(field[1, :, :]),
            "time" => "NP1",
        )
    end

end

function switch_bonds!(field_N, field_NP1)
    for fieldID in eachindex(field_NP1)
        copyto!(field_N[fieldID], field_NP1[fieldID])
    end
end

"""
    switch_NP1_to_N()

Switches the fields from NP1 to N.active
"""
function switch_NP1_to_N()
    active = _get_field("Active")
    for key in keys(data["NP1_to_N"])
        if key == "Bond Damage"
            field_NP1 = get_field(key, "NP1")
            field_N = get_field(key, "N")
            switch_bonds!(field_N, field_NP1)
            continue
        end
        data["NP1_to_N"][key][1], data["NP1_to_N"][key][2] =
            data["NP1_to_N"][key][2], data["NP1_to_N"][key][1]
        field_NP1 = get_field(key, "NP1")
        if field_NP1[1] isa AbstractVector || field_NP1[1] isa AbstractArray
            fill_in_place!(field_NP1, data["NP1_to_N"][key][3], active)
        else
            fill!(field_NP1, data["NP1_to_N"][key][3])
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
    synch_fields = get_synch_fields()
    # @debug synch_fields
    for synch_field in keys(synch_fields)
        synchronise_field(
            get_comm(),
            synch_fields,
            get_overlap_map(),
            get_field,
            synch_field,
            direction,
        )
    end
end
end
