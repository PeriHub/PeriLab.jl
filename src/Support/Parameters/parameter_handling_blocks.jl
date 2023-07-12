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


