function get_horizon(data, blockID)
    check = check_element(data["Blocks"], "block_" * string(blockID))
    if check
        check2 = check_element(data["Blocks"]["block_"*string(blockID)], "Horizon")
        if check2
            return data["Blocks"]["block_"*string(blockID)]["Horizon"]
        else
            @error "Horizon of Block $blockID is not defined"
        end
        @error "Block $blockID is not defined"
    end
end

function get_number_of_blocks(data)
    check = check_element(data, "Blocks")
    if check
        return length(data["Blocks"])
    else
        @error "Block $blockID is not defined"
    end
    return 0
end