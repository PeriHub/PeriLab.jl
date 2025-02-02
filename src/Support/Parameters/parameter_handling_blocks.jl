# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export get_density
export get_fem_block
export get_heat_capacity
export get_horizon
export get_values
export get_number_of_blocks
export get_block_models
export get_angles

"""
    get_block_names(params::Dict, block_ids::Vector{Int64})

Get the names of the blocks.

# Arguments
- `params::Dict`: The parameters dictionary.
- `block_ids::Vector{Int64}`: The IDs of the blocks
# Returns
- `block_names::Vector{String}`: The names of the blocks.
"""
function get_block_names(params::Dict, block_ids::Vector{Int64})
    param_block_ids = [v["Block ID"] for v in values(params["Blocks"])]
    block_list = Vector{String}()
    for id = 1:maximum(param_block_ids)
        if !(id in block_ids)
            @warn "Block with ID $id is not defined in the provided mesh"
            continue
        end
        for (block, value) in params["Blocks"]
            if value["Block ID"] == id
                push!(block_list, block)
            end
        end
    end
    return block_list
end

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
    get_fem_block(params::Dict, block_id::Int64)

Get the fem_block of a block.

# Arguments
- `params::Dict`: The parameters
- `block_id::Int64`: The ID of the block
# Returns
- `fem_block::Float64`: The fem_block of the block
"""
function get_fem_block(params::Dict, block_id::Int64)
    return get_values(params, block_id, "FEM", false)
end

"""
    get_heat_capacity(params::Dict, block_id::Int64)

Get the heat capacity of a block.

# Arguments
- `params::Dict`: The parameters
- `block_id::Int64`: The ID of the block
# Returns
- `heat_capacity::Float64`: The heat capacity of the block
"""
function get_heat_capacity(params::Dict, block_id::Int64)
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
    get_angles(params::Dict, block_id::Int64, dof::Int64)

Get the horizon of a block.

# Arguments
- `params::Dict`: The parameters
- `block_id::Int64`: The ID of the block
- `dof::Int64`: The dof
# Returns
- `angles::Float64`: The angles of the block
"""
function get_angles(params::Dict, block_id::Int64, dof::Int64)

    if !haskey(params["Blocks"]["block_"*string(block_id)], "Angle X")
        return nothing
    end
    if dof == 2
        return get_values(params, block_id, "Angle X")
    elseif dof == 3
        return [
            get_values(params, block_id, "Angle X"),
            get_values(params, block_id, "Angle Y"),
            get_values(params, block_id, "Angle Z"),
        ]
    end
end

"""
    get_values(params::Dict, block_name::String, valueName::String, defaultValue::Union{Float64,Bool,Nothing})

Get the value of a block.

# Arguments
- `params::Dict`: The parameters
- `block_name::String`: The ID of the block
- `valueName::String`: The name of the value
- `defaultValue::Union{Float64,Bool,Nothing`: The default value
# Returns
- `value::Float64`: The value of the block
"""
function get_values(
    params::Dict,
    block_name::String,
    valueName::String,
    defaultValue::Union{Float64,Bool,Nothing} = nothing,
)
    if haskey(params["Blocks"], block_name)
        if haskey(params["Blocks"][block_name], valueName)
            return params["Blocks"][block_name][valueName]
        end
        if isnothing(defaultValue)
            @error "$valueName of $block_name is not defined"
        end
        return defaultValue
    end
    @error "$block_name is not defined"
    return nothing
end

function get_values(
    params::Dict,
    block_id::Int64,
    valueName::String,
    defaultValue::Union{Float64,Bool,Nothing} = nothing,
)
    for block_name in keys(params["Blocks"])
        if params["Blocks"][block_name]["Block ID"] == block_id
            if haskey(params["Blocks"][block_name], valueName)
                return params["Blocks"][block_name][valueName]
            end
            if isnothing(defaultValue)
                @error "$valueName of $block_name is not defined"
            end
            return defaultValue
        end
    end
    @error "$block_name is not defined"
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
