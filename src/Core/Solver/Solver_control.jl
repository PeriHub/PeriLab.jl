
include("../../Support/Parameters/parameter_handling.jl")
include("../../Physics/Physics_Factory.jl")
import .Physics
module Solver


function init(params, datamodel)

    blockNodes = get_blockNodes(datamodel.get_field("Block_Id"))
    mechanical, thermal, additive = get_solver_optopms(params)
    if mechanical
        dof = datamanager.get_dof()
        force = datamanager.create_node_field("Force", Float32, dof)
        Y = datamanager.create_node_field("Deformed State", Float32, dof)
        u = datamanager.create_node_field("Displacement", Float32, dof)
        bu = datamanager.create_bond_field("Deformed Bond Geometry", Float32, dof + 1)
    end
    if thermal
        temperature = datamanager.create_node_field("Temperature", Float32, 1)
        flow = datamanager.create_node_field("Flow", Float32, dof) # -> check dof
    end
    if additive

    end
    physics = Physics.get_physics(params)


    return blockNodes, datamodel
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