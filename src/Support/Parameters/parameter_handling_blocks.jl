# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    get_density(params::Dict, block_id::Int64)

    Get the density of a block.

    # Arguments
   - `params::Dict`: The parameters
   - `block_id::Int64`: The ID of the block
   # Returns
   - `density::Float64`: The density of the block
"""
function get_density(params::Dict, block_id::Int64)
    return get_values(params, block_id, "Density")
end

"""
    get_heatcapacity(params::Dict, block_id::Int64)

    Get the heat capacity of a block.

    # Arguments
   - `params::Dict`: The parameters
   - `block_id::Int64`: The ID of the block
   # Returns
   - `heatcapacity::Float64`: The heat capacity of the block
"""
function get_heatcapacity(params::Dict, block_id::Int64)
    return get_values(params, block_id, "Specific Heat Capacity")
end

"""
    get_horizon(params::Dict, block_id::Int64)

    Get the horizon of a block.

    # Arguments
   - `params::Dict`: The parameters
   - `block_id::Int64`: The ID of the block
   # Returns
   - `horizon::Float64`: The horizon of the block
"""
function get_horizon(params::Dict, block_id::Int64)
    return get_values(params, block_id, "Horizon")
end

"""
    get_values(params::Dict, block_id::Int64, valueName::String)

    Get the value of a block.

    # Arguments
    - `params::Dict`: The parameters
    - `block_id::Int64`: The ID of the block
    - `valueName::String`: The name of the value
    # Returns
    - `value::Float64`: The value of the block
"""
function get_values(params::Dict, block_id::Int64, valueName::String)
    if haskey(params["Blocks"], "block_" * string(block_id))
        if haskey(params["Blocks"]["block_"*string(block_id)], valueName)
            return params["Blocks"]["block_"*string(block_id)][valueName]
        end
        @error "$valueName of Block $block_id is not defined"
        return nothing
    end
    @error "Block $block_id is not defined"
    return nothing
end

"""
    get_number_of_blocks(params::Dict)

    Get the number of blocks.

    # Arguments
    - `params::Dict`: The parameters
    # Returns
    - `number_of_blocks::Int64`: The number of blocks
"""
function get_number_of_blocks(params::Dict)
    check = haskey(params::Dict, "Blocks")
    if haskey(params::Dict, "Blocks") && length(params["Blocks"]) > 0
        return length(params["Blocks"])
    end
    @error "No blocks defined"
    return nothing
end

"""
    get_block_models(params::Dict, block_id::Int64)

    Get the models of a block.

    # Arguments
    - `params::Dict`: The parameters
    - `block_id::Int64`: The ID of the block
    # Returns
    - `modelDict::Dict{String,String}`: The models of the block
"""
function get_block_models(params::Dict, block_id::Int64)
    modelDict = Dict{String,String}()
    check = haskey(params["Blocks"], "block_" * string(block_id))
    if check
        for key in keys(params["Blocks"]["block_"*string(block_id)])
            if occursin("Model", key)
                modelDict[key] = params["Blocks"]["block_"*string(block_id)][key]
            end
        end
    end
    return modelDict
end

