# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function sum_block(datamanager::Module, fieldKey::String, block::Int64, blockNodes::Dict{Int64,Vector{Int64}})
    # get blockNodes
    # check NP1
    if fieldKey * "NP1" in datamanager.get_all_field_keys()
        fieldKey *= "NP1"
    end
    if !(fieldKey in datamanager.get_all_field_keys())
        @error "Field $fieldKey does not exists for compute sum."
        return Nothing
    end
    field = datamanager.get_field(fieldKey)
    nodes = blockNodes[block]

    return global_value_sum(field, nodes)
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
