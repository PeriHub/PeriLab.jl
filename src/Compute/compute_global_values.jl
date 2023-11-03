# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function calculate_nodelist(datamanager::Module, fieldKey::String, calculationType::String, node_set::Vector{Int64})
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

    if calculationType == "Sum"
        return global_value_sum(field, nodes)
    elseif calculationType == "Maximum"
        return global_value_max(field, nodes)
    elseif calculationType == "Minimum"
        return global_value_min(field, nodes)
    elseif calculationType == "Average"
        return global_value_avg(field, nodes)
    else
        @warn "Unknown calculation type $calculationType"
        return nothing
    end
end

function calculate_block(datamanager::Module, fieldKey::String, calculationType::String, block::Int64)
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
    blockIDs = datamanager.get_field("Block_Id")
    nnodes = datamanager.get_nnodes()
    blockNodes = findall(item -> item == block, blockIDs[1:nnodes])

    if calculationType == "Sum"
        return global_value_sum(field, blockNodes)
    elseif calculationType == "Maximum"
        return global_value_max(field, blockNodes)
    elseif calculationType == "Minimum"
        return global_value_min(field, blockNodes)
    elseif calculationType == "Average"
        return global_value_avg(field, blockNodes)
    else
        @warn "Unknown calculation type $calculationType"
        return 0
    end
end

function global_value_sum(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = sum(field[nodes, iID])
    end
    return returnValue
end

function global_value_max(field::SubArray, nodes::Union{SubArray,Vector{Int64}})
    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = maximum(field[nodes, iID])
    end
    return returnValue
end

function global_value_min(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = minimum(field[nodes, iID])
    end
    return returnValue
end

function global_value_avg(field::SubArray, nodes::Union{SubArray,Vector{Int64}})

    returnValue = zeros(length(field[1, :]))
    for iID in eachindex(field[1, :])
        returnValue[iID] = sum(field[nodes, iID]) / length(nodes)
    end
    return returnValue
end
