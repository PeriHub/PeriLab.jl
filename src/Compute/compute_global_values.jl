# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    calculate_nodelist(datamanager::Module, fieldKey::String, calculation_type::String, node_set::Vector{Int64})

    Calculate the global value of a field for a given set of nodes.

    # Arguments
   - `datamanager::Data_manager`: Datamanager.
   - `fieldKey::String`: Field key.
   - `calculation_type::String`: Calculation type.
   - `node_set::Vector{Int64}`: Node set.
    # Returns
   - `value::Vector`: Global value.
   - `nnodes::Int64`: Number of nodes.
"""
function calculate_nodelist(datamanager::Module, fieldKey::String, calculation_type::String, node_set::Vector{Int64})
    # get blockNodes
    # check NP1
    if fieldKey * "NP1" in datamanager.get_all_field_keys()
        fieldKey *= "NP1"
    end
    if !(fieldKey in datamanager.get_all_field_keys())
        @error "Field $fieldKey does not exists for compute sum."
        return nothing
    end
    field = datamanager.get_field(fieldKey)
    nnodes = datamanager.get_nnodes()
    nodes = 1:nnodes
    nodes = nodes[node_set]
    if calculation_type == "Sum"
        if length(nodes) == 0
            value = fill(field_type(0), length(field[1, :]))
        end
        value = global_value_sum(field, nodes)
    elseif calculation_type == "Maximum"
        if length(nodes) == 0
            value = fill(typemin(field_type), length(field[1, :]))
        end
        value = global_value_max(field, nodes)
    elseif calculation_type == "Minimum"
        if length(nodes) == 0
            value = fill(typemax(field_type), length(field[1, :]))
        end
        value = global_value_min(field, nodes)
    elseif calculation_type == "Average"
        if length(nodes) == 0
            value = fill(field_type(0), length(field[1, :]))
        end
        value = global_value_avg(field, nodes)
    else
        @warn "Unknown calculation type $calculationType"
        return nothing
    end
    return value, Int64(length(nodes))
end

"""
    calculate_block(datamanager::Module, fieldKey::String, calculation_type::String, block::Int64)

    Calculate the global value of a field for a given block.

    # Arguments
   - `datamanager::Data_manager`: Datamanager.
   - `fieldKey::String`: Field key.
   - `calculation_type::String`: Calculation type.
   - `block::Int64`: Block number.
    # Returns
   - `value::Vector`: Global value.
   - `nnodes::Int64`: Number of nodes.
"""
function calculate_block(datamanager::Module, fieldKey::String, calculation_type::String, block::Int64)
    # get blockNodes
    # check NP1
    if fieldKey * "NP1" in datamanager.get_all_field_keys()
        fieldKey *= "NP1"
    end
    if !(fieldKey in datamanager.get_all_field_keys())
        @error "Field $fieldKey does not exists for compute sum."
        return nothing
    end
    field = datamanager.get_field(fieldKey)
    field_type = datamanager.get_field_type(fieldKey)
    block_ids = datamanager.get_field("Block_Id")
    nnodes = datamanager.get_nnodes()
    blockNodes = findall(item -> item == block, block_ids[1:nnodes])

    if calculation_type == "Sum"
        if length(blockNodes) == 0
            value = fill(field_type(0), length(field[1, :]))
        else
            value = global_value_sum(field, blockNodes)
        end
    elseif calculation_type == "Maximum"
        if length(blockNodes) == 0
            value = fill(typemin(field_type), length(field[1, :]))
        else
            value = global_value_max(field, blockNodes)
        end
    elseif calculation_type == "Minimum"
        if length(blockNodes) == 0
            value = fill(typemax(field_type), length(field[1, :]))
        else
            value = global_value_min(field, blockNodes)
        end
    elseif calculation_type == "Average"
        if length(blockNodes) == 0
            value = fill(field_type(0), length(field[1, :]))
        else
            value = global_value_avg(field, blockNodes)
        end
    else
        @warn "Unknown calculation type $calculation_type"
        return nothing
    end
    return value, Int64(length(blockNodes))
end

"""
    global_value_sum(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    Calculate the global sum of a field for given nodes.
    
    # Arguments
   - `field::SubArray`: Field.
   - `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
    # Returns
   - `returnValue::Vector`: Global value.
"""
function global_value_sum(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = sum(field[nodes, iID])
    end
    return returnValue
end

"""
    global_value_max(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    Calculate the global maximum of a field for given nodes.
    
    # Arguments
   - `field::SubArray`: Field.
   - `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
    # Returns
   - `returnValue::Vector`: Global value.
"""
function global_value_max(field::SubArray, nodes::Union{SubArray,Vector{Int64}})
    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = maximum(field[nodes, iID])
    end
    return returnValue
end

"""
    global_value_min(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    Calculate the global minimum of a field for given nodes.

    # Arguments
   - `field::SubArray`: Field.
   - `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
    # Returns
   - `returnValue::Vector`: Global value.
"""
function global_value_min(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = minimum(field[nodes, iID])
    end
    return returnValue
end

"""
    global_value_avg(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    Calculate the global average of a field for given nodes.
    
    # Arguments
   - `field::SubArray`: Field.
   - `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
    # Returns
   - `returnValue::Vector`: Global value.
"""
function global_value_avg(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = sum(field[nodes, iID]) / length(nodes)
    end
    return returnValue
end
