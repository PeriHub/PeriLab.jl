# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_density(params::Dict, blockID::Int64)
    return get_values(params, blockID, "Density")
end

function get_heatcapacity(params::Dict, blockID::Int64)
    return get_values(params, blockID, "Specific Heat Capacity")
end

function get_horizon(params::Dict, blockID::Int64)
    return get_values(params, blockID, "Horizon")
end

function get_values(params::Dict, blockID::Int64, valueName::String)
    if check_element(params["Blocks"], "block_" * string(blockID))
        if check_element(params["Blocks"]["block_"*string(blockID)], valueName)
            return params["Blocks"]["block_"*string(blockID)][valueName]
        end
        @error "$valueName of Block $blockID is not defined"
        return
    end
    @error "Block $blockID is not defined"
    return
end

function get_number_of_blocks(params::Dict)
    check = check_element(params::Dict, "Blocks")
    if check_element(params::Dict, "Blocks") && length(params["Blocks"]) > 0
        return length(params["Blocks"])
    end
    @error "No blocks defined"
    return
end

function get_block_models(params::Dict, blockID::Int64)
    modelDict = Dict{String,String}()
    check = check_element(params["Blocks"], "block_" * string(blockID))
    if check
        for key in keys(params["Blocks"]["block_"*string(blockID)])
            if occursin("Model", key)
                modelDict[key] = params["Blocks"]["block_"*string(blockID)][key]
            end
        end
    end
    return modelDict
end

