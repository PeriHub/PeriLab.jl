
include("../../Support/Parameters/parameter_handling.jl")
include("../../Physics/Physics_Factory.jl")
import .Physics
module Solver


function init(params, datamanager)

    blockNodes = get_blockNodes(datamodel.get_field("Block_Id"))
    mechanical, thermal, additive = get_solver_options(params)
    if mechanical
        dof = datamanager.get_dof()
        force = datamanager.create_node_field("Force", Float32, dof)
        Y = datamanager.create_node_field("Deformed State", Float32, dof)
        u = datamanager.create_node_field("Displacements", Float32, dof)
        bu = datamanager.create_bond_field("Deformed Bond Geometry", Float32, dof + 1)
    end
    if thermal
        temperature = datamanager.create_node_field("Temperature", Float32, 1)
        flow = datamanager.create_node_field("Flow", Float32, dof) # -> check dof
    end
    if additive

    end
    physics = Physics.get_physics(params)
    boundary_condition(params, datamanager)

    return blockNodes, datamanager
end

function boundary_condition(params, datamanager)
    global_nset = get_node_sets(params)
    bcs = et_node_bcs(params)
    return datamanager.glob_to_loc(global_nset), bcs
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

function run_solver()

end

function solver()
    blockNodes, datamodel = init(params, datamodel)
    run_solver(blockNodes, datamodel)
end

end