# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Data_Manager
using MPI
using DataStructures: OrderedDict

##########################
# Variables
##########################

const fields = Dict()
const data = Dict()
#####################

struct DataField{T,N}
    name::String
    data::Array{T,N}
    bond_or_node::String
end

struct FieldManager
    fields::Dict{String,DataField}
end

mutable struct NP1_to_N{T}
    N::String
    NP1::String
    value::T
end

fieldmanager = FieldManager(Dict{String,DataField}())
#####

include("./Data_manager/data_manager_contact.jl")
include("./Data_manager/data_manager_FEM.jl")
include("./Data_manager/data_manager_fields.jl")
include("./Data_manager/data_manager_MPI.jl")
include("./Data_manager/data_manager_solving_process.jl")
include("./Data_manager/data_manager_status_vars.jl")
include("./Data_manager/data_manager_utilities.jl")
export add_active_model
export fem_active
export initialize_data
export get_active_models
export get_all_field_keys
export get_bc_free_dof
export get_accuracy_order
export get_aniso_crit_values
export get_block_name_list
export get_block_id_list
export get_crit_values_matrix
export get_comm
export get_coupling_dict
export get_coupling_fe_nodes
export get_directory
export get_field
export get_field_type
export get_inverse_nlist
export get_local_nodes
export get_model_module
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
export get_output_frequency
export get_element_rotation
export init_properties
export remove_active_model
export set_accuracy_order
export set_bc_free_dof
export set_block_name_list
export set_block_id_list
export set_crit_values_matrix
export set_coupling_dict
export set_coupling_fe_nodes
export set_aniso_crit_values
export set_directory
export set_inverse_nlist
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

"""
    initialize_data()

Initialize all parameter in the Data_Manager and sets them to the default values.
"""
function initialize_data()
    data["current_time"] = 0.0

    data["step"] = 0
    data["max_step"] = 0
    data["nnodes"] = 0
    data["num_controller"] = 0
    data["num_responder"] = 0
    data["num_elements"] = 0
    data["nnsets"] = 0
    data["dof"] = 0
    data["fem_option"] = false
    data["block_name_list"] = Vector{String}()
    data["block_id_list"] = Vector{Int64}()
    data["distribution"] = []
    data["crit_values_matrix"] = fill(-1, (1, 1, 1))
    data["aniso_crit_values"] = Dict()
    data["properties"] = OrderedDict()
    data["glob_to_loc"] = Dict()
    #data["field_array_type"] = Dict()
    data["field_types"] = Dict()
    data["field_names"] = Vector{String}([])
    data["fields_to_synch"] = Dict()
    data["local_fields_to_synch"] = Dict("Material Model" => Dict(),
                                         "Damage Model" => Dict(),
                                         "Thermal Model" => Dict(),
                                         "Pre Calculation Model" => Dict(),
                                         "Degradation Model" => Dict(),
                                         "Additive Model" => Dict(),
                                         "Surface Correction" => Dict())
    data["filedirectory"] = ""
    data["inverse_nlist"] = []
    data["model_modules"] = OrderedDict{String,Module}()
    data["analysis_models"] = Dict()
    data["nsets"] = Dict{String,Vector{Int64}}()
    data["overlap_map"] = Dict{Int64,Dict{Int64,Dict{String,Vector{Int64}}}}()
    data["contact_overlap_map"] = Dict()
    data["pre_calculation_order"] = [
        "Deformed Bond Geometry",
        "Shape Tensor",
        "Deformation Gradient",
        "Bond Associated Correspondence"
    ]
    data["coupling_dict"] = Dict{Int64,Int64}()
    data["output_frequency"] = []
    data["accuracy_order"] = 1
    data["rank"] = 0
    data["cancel"] = false
    data["max_rank"] = 0
    data["mpi_active"] = false
    data["silent"] = false
    data["verbose"] = false
    data["rotation"] = false
    data["element_rotation"] = false
    data["active_models"] = OrderedDict{String,Module}()
    data["all_active_models"] = OrderedDict{String,Module}()
    data["NP1_to_N"] = Dict{String,NP1_to_N}()
    data["coupling_fe_nodes"] = []
    data["BC_free_dof"] = []
    data["Contact Nodes"] = Dict{Int64,Vector{Int64}}()
    data["All positions"] = []
    data["All Blocks"] = []
    data["Free Surfaces"] = Dict()
    data["Free Surface Connections"] = Dict{Int64,Vector{Int64}}() # point ID -> surface IDs
    data["Free Surface Nodes"] = Dict{Int64,Vector{Int64}}() # block ID -> surface point IDs
    data["Contact Dictionary"] = Dict()
    data["Global Contact IDs"] = Vector{Int64}([])
    data["Local Contact IDs"] = Dict{Int64,Int64}()
    data["Contact Properties"] = Dict()
    data["Contact block IDs"] = Dict()
    data["Exchange id to local id"] = Dict{Int64,Int64}()
    data["Contact Search Step"] = Dict{Any,Int64}()
    data["Contact Search No Pairs"] = Dict{Any,Bool}()
    data["Synchronization List"] = Dict{Int64,Int64}()
    data["Global Master Search Nodes"] = Dict()
    data["Global Slave Search Nodes"] = Dict()
    fields[Int64] = Dict()
    fields[Float64] = Dict()
    fields[Bool] = Dict()
end
###################################

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
    get_block_name_list()

Returns a list of all block IDs.
"""
function get_block_name_list()
    return data["block_name_list"]
end

"""
    get_block_id_list()

Returns a list of all block IDs.
"""
function get_block_id_list()
    return data["block_id_list"]
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
    get_field_type()
Get the type of a field

# Returns
- `get_field_type` (string): returns the type of a field
"""
function get_field_type(name::String, vartype::Bool = true)
    if name in data["field_names"]
        if vartype
            return data["field_types"][name]["vartype"]
        else
            return data["field_types"][name]["type"]
        end
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
    return [data["glob_to_loc"][global_node]
            for
            global_node in global_nodes
            if global_node in keys(data["glob_to_loc"])]
end

function get_model_module(entry::AbstractString)
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
function get_properties(block_id::Int64, property::String)::Dict{String,Any}
    if check_property(block_id, property)
        return convert(Dict{String,Any}, data["properties"][block_id][property]) # TODO check why it is needed!
    end
    return Dict{String,Any}()::Dict{String,Any}
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
function get_property(block_id::Int64, property::String,
                      value_name::String)
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
function get_rank()::Int64
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
function get_max_rank()::Int64
    return data["max_rank"]
end

"""
    get_mpi_active()

This function returns if MPI is active.

# Returns
- `mpi_active`::Bool

# Example
```julia
rank = get_mpi_active()
"""
function get_mpi_active()::Bool
    return data["mpi_active"]
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
    get_verbose()

This function returns the `verbose` flag.

# Returns
- `verbose`::Bool: The value of the `verbose` variable.
"""
function get_verbose()
    return data["verbose"]
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
    loc_to_glob(range::UnitRange{Int64})

Converts the local index to the global index.

# Arguments
- `range::Union{UnitRange{Int64}, Vector{Int64}}`: The range of the local index.

Example:
```julia
loc_to_glob(1:10)  # converts the local index to the global index
```
"""
function loc_to_glob(range::Vector{Int64})
    return data["distribution"][range]
end
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
    block_id_list = get_block_id_list()
    for iblock in block_id_list
        data["properties"][iblock] = OrderedDict{String,Dict{String,Any}}()

        for prop_name in ["Additive Model", "Damage Model", "Pre Calculation Model",
            "Thermal Model", "Degradation Model", "Material Model"]
            data["properties"][iblock][prop_name] = Dict{String,Any}()
        end
    end
    return collect(keys(data["properties"][block_id_list[1]]))
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
    set_block_name_list(blocks_name_list::Vector{String})

Sets the block list globally.

# Arguments
- `blocks_name_list::Vector{String}`: The block list.
"""
function set_block_name_list(blocks_name_list::Vector{String})
    data["block_name_list"] = blocks_name_list
end

"""
    set_block_id_list(blocks_id_list::Vector{Int64})

Sets the block list globally.

# Arguments
- `blocks_id_list::Vector{Int64}`: The block list.
"""
function set_block_id_list(blocks_id_list::Vector{Int64})
    data["block_id_list"] = blocks_id_list
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
    set_inverse_nlist(inv_nlist::Vector{Dict{Int64,Int64}})

Sets the inverse nlist globally.

# Arguments
- `inv_nlist::Vector{Dict{Int64,Int64}}`: The inverse nlist.
"""
function set_inverse_nlist(inv_nlist::Vector{Dict{Int64,Int64}})
    data["inverse_nlist"] = inv_nlist
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
function set_property(block_id::Int64, property::String, value_name::String,
                      value::T) where {T}
    prop_dict = get!(data["properties"][block_id], property, Dict{String,Any}())
    prop_dict[value_name] = value
end

function set_property(property::String, value_name::String, value::T) where {T}
    for block_id in eachindex(data["properties"])
        set_property(block_id, property, value_name, value)
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
function set_properties(block_id::Int64, property::String, values::T) where {T}
    if values isa Dict{String,Any}
        data["properties"][block_id][property] = values
    else
        data["properties"][block_id][property] = convert(Dict{String,Any}, values)
    end
end

"""
    set_properties(property, values)

Sets the values of a specified `property` for a all `blocks`. E.g. for FEM, because it corresponds not to a block yet,

# Arguments
- `property`::String: The name of the property.
- `values`::Any: The values to set for the specified `property`.
"""
function set_properties(property::String, values::T) where {T}
    for id in eachindex(data["properties"])
        # Typ-sichere Zuweisung
        if values isa Dict{String,Any}
            data["properties"][id][property] = values
        else
            data["properties"][id][property] = convert(Dict{String,Any}, values)
        end
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

function set_model_module(entry::AbstractString, mod::Module)
    data["model_modules"][entry] = mod
end

"""
    set_max_rank(value::Int64)

Sets the maximum rank globally.

# Arguments
- `value::Int64`: The value to set as the maximum rank.
"""
function set_max_rank(value::Int64)
    if value > 1
        data["mpi_active"] = true
    end
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
    set_silent(value::Bool)

Sets the silent flag.

# Arguments
- `value::Bool`: The silent flag.
"""
function set_silent(value::Bool)
    data["silent"] = value
end

"""
    set_verbose(value::Bool)

Sets the verbose flag.

# Arguments
- `value::Bool`: The verbose flag.
"""
function set_verbose(value::Bool)
    data["verbose"] = value
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
        data["fields_to_synch"][name] = Dict{String,Any}("upload_to_cores" => upload_to_cores,
                                                         "download_from_cores" => download_from_cores,
                                                         "dof" => dof,
                                                         "time" => "Constant")
    elseif name * "NP1" in get_all_field_keys()
        field = get_field(name, "NP1")
        if dof == 0
            dof = length(field[1, :, :])
        end
        data["fields_to_synch"][name] = Dict{String,Any}("upload_to_cores" => upload_to_cores,
                                                         "download_from_cores" => download_from_cores,
                                                         "dof" => length(field[1, :, :]),
                                                         "time" => "NP1")
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
        @error "Model $model is not defined. If it is a new model type, please add it to data[\"local_fields_to_synch\"] in the Data_Manager."
        return nothing
    end
    if name in get_all_field_keys()
        field = _get_field(name)
        if dof == 0
            dof = length(field[1, :, :])
        end
        data["local_fields_to_synch"][model][name] = Dict{String,Any}("upload_to_cores" => upload_to_cores,
                                                                      "download_from_cores" => download_from_cores,
                                                                      "dof" => dof,
                                                                      "time" => "Constant")
    elseif name * "NP1" in get_all_field_keys()
        field = get_field(name, "NP1")
        if dof == 0
            dof = length(field[1, :, :])
        end
        data["local_fields_to_synch"][model][name] = Dict{String,Any}("upload_to_cores" => upload_to_cores,
                                                                      "download_from_cores" => download_from_cores,
                                                                      "dof" => length(field[1,
                                                                                            :,
                                                                                            :]),
                                                                      "time" => "NP1")
    end
end
end
