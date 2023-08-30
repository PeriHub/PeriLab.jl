
function get_density(params, blockID)
    if check_element(params["Blocks"], "block_" * string(blockID))
        return params["Blocks"]["block_"*string(blockID)]["Density"]
    end
    @error "No density defined for block " * string(blockID)
    return
end

function get_horizon(params, blockID)
    check = check_element(params["Blocks"], "block_" * string(blockID))
    if check
        check2 = check_element(params["Blocks"]["block_"*string(blockID)], "Horizon")
        if check2
            return params["Blocks"]["block_"*string(blockID)]["Horizon"]
        else
            @error "Horizon of Block $blockID is not defined"
        end
        @error "Block $blockID is not defined"
    end
end

function get_number_of_blocks(params)
    check = check_element(params, "Blocks")
    if check
        return length(params["Blocks"])
    else
        @error "No blocks defined"
    end
    return 0
end

function get_block_models(params, blockID)
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

