
include("../../Physics/Physics_Factory.jl")
import .Physics
module Solver


function init(params, datamodel)

    blockNodes = get_blockNodes(datamodel.get_field("Block_Id"))
    physics = Physics.get_physics(params)


    return blockNodes
end
function get_blockNodes(blockID)
    maxBlock = maximum(blockID)
    blockNodes = distribution = [collect(1:maxBlock)]
    for i in 1:maxBlock
        blockNodes[i] = find_indices(blockID, i)
    end
    return blockNodes
end


function solver(params, datamodel)
    blockNodes = init(params, datamodel)
end


end