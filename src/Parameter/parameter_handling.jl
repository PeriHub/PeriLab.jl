using MPI
function get_mesh_name(data)
    check = check_element(data["Discretization"], "Input Mesh File")
    if !check
        @error "No mesh file is defined."
    end
    return data["Discretization"]["Input Mesh File"]
end

function get_topology_name(data)
    check = check_element(data["Discretization"], "Input FEM Topology File")
    return check, data["Discretization"]["Input FEM Topology File"]
end

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


function check_element(data, key)
    return haskey(data, key)
end