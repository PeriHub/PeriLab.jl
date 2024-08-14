# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    calculate_nodelist(datamanager::Module, field_key::String, dof::Union{Int64,Vector{Int64}}, calculation_type::String, node_set::Vector{Int64})

Calculate the global value of a field for a given set of nodes.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `field_key::String`: Field key.
- `dof::Union{Int64,Vector{Int64}}`: Degree of freedom
- `calculation_type::String`: Calculation type.
- `node_set::Vector{Int64}`: Node set.
# Returns
- `value::Vector`: Global value.
- `nnodes::Int64`: Number of nodes.
"""
function calculate_nodelist(
    datamanager::Module,
    field_key::String,
    dof::Union{Int64,Vector{Int64}},
    calculation_type::String,
    node_set::Vector{Int64},
)
    # get block_nodes
    # check NP1
    if datamanager.has_key(field_key * "NP1")
        field_key *= "NP1"
    end
    if !datamanager.has_key(field_key)
        @error "Field $field_key does not exists for compute sum."
        return nothing
    end
    field = datamanager.get_field(field_key)
    field_type = datamanager.get_field_type(field_key)
    node_list = datamanager.get_local_nodes(node_set)

    if calculation_type == "Sum"
        if length(node_list) == 0
            value = field_type(0)
        else
            value = global_value_sum(field, dof, node_list)
        end
    elseif calculation_type == "Maximum"
        if length(node_list) == 0
            value = typemin(field_type)
        else
            value = global_value_max(field, dof, node_list)
        end
    elseif calculation_type == "Minimum"
        if length(node_list) == 0
            value = typemax(field_type)
        else
            value = global_value_min(field, dof, node_list)
        end
    elseif calculation_type == "Average"
        if length(node_list) == 0
            value = field_type(0)
        else
            value = global_value_avg(field, dof, node_list)
        end
    else
        @warn "Unknown calculation type $calculation_type"
        return nothing
    end
    return value, Int64(length(node_list))
end

"""
    calculate_block(datamanager::Module, field_key::String, dof::Int64, calculation_type::String, block::Int64)

Calculate the global value of a field for a given block.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `field_key::String`: Field key.
- `dof::Union{Int64,Vector{Int64}}`: Degree of freedom
- `calculation_type::String`: Calculation type.
- `block::Int64`: Block number.
# Returns
- `value::Float64`: Global value.
- `nnodes::Int64`: Number of nodes.
"""
function calculate_block(
    datamanager::Module,
    field_key::String,
    dof::Union{Int64,Vector{Int64}},
    calculation_type::String,
    block::Int64,
)
    # get block_nodes
    block_ids = datamanager.get_field("Block_Id")
    nnodes = datamanager.get_nnodes()
    block_nodes = Vector{Int64}(findall(item -> item == block, block_ids[1:nnodes]))
    return calculate_nodelist(datamanager, field_key, dof, calculation_type, block_nodes)
end

"""
    global_value_sum(field::SubArray, dof::Union{Int64,Vector{Int64}}, nodes::Union{SubArray,Vector{Int64}})

Calculate the global sum of a field for given nodes.

# Arguments
- `field::SubArray`: Field.
- `dof::Union{Int64,Vector{Int64}`: Degree of freedom
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
# Returns
- `returnValue::Vector`: Global value.
"""
function global_value_sum(
    field::SubArray,
    dof::Union{Int64,Vector{Int64}},
    nodes::Union{SubArray,Vector{Int64}},
)

    if dof isa Int64
        return sum(field[nodes, dof])
    else
        return sum(field[nodes, dof[1], dof[2]])
    end
end

"""
    global_value_max(field::SubArray, dof::Union{Int64,Vector{Int64}}, nodes::Union{SubArray,Vector{Int64}})

Calculate the global maximum of a field for given nodes.

# Arguments
- `field::SubArray`: Field.
- `dof::Union{Int64,Vector{Int64}}`: Degree of freedom
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
# Returns
- `returnValue::Vector`: Global value.
"""
function global_value_max(
    field::SubArray,
    dof::Union{Int64,Vector{Int64}},
    nodes::Union{SubArray,Vector{Int64}},
)
    # returnValue = zeros(length(field[1, :]))
    # for iID in eachindex(field[1, :])
    #     returnValue[iID] = maximum(field[nodes, iID])
    # end
    # return returnValue
    return maximum(field[nodes, dof])
end

"""
    global_value_min(field::SubArray, dof::Union{Int64,Vector{Int64}}, nodes::Union{SubArray,Vector{Int64}})

Calculate the global minimum of a field for given nodes.

# Arguments
- `field::SubArray`: Field.
- `dof::Union{Int64,Vector{Int64}}`: Degree of freedom
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
# Returns
- `returnValue::Vector`: Global value.
"""
function global_value_min(
    field::SubArray,
    dof::Union{Int64,Vector{Int64}},
    nodes::Union{SubArray,Vector{Int64}},
)

    # returnValue = zeros(length(field[1, :]))
    # for iID in eachindex(field[1, :])
    #     returnValue[iID] = minimum(field[nodes, iID])
    # end
    # return returnValue
    return minimum(field[nodes, dof])
end

"""
    global_value_avg(field::SubArray, dof::Union{Int64,Vector{Int64}}, nodes::Union{SubArray,Vector{Int64}})

Calculate the global average of a field for given nodes.

# Arguments
- `field::SubArray`: Field.
- `dof::Union{Int64,Vector{Int64}}`: Degree of freedom
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
# Returns
- `returnValue::Vector`: Global value.
"""
function global_value_avg(
    field::SubArray,
    dof::Union{Int64,Vector{Int64}},
    nodes::Union{SubArray,Vector{Int64}},
)

    # returnValue = zeros(length(field[1, :]))
    # for iID in eachindex(field[1, :])
    #     returnValue[iID] = sum(field[nodes, iID]) / length(nodes)
    # end
    # return returnValue
    return sum(field[nodes, dof]) / length(nodes)
end
