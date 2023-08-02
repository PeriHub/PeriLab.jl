
include("../../Support/Parameters/parameter_handling.jl")
include("../../Physics/Physics_Factory.jl")
import .Physics
module Solver


function init(params, datamanager)
    output_filenames = get_output_filenames(params)

    init_write_results(output_filenames, datamanager)
    blockNodes = get_blockNodes(datamanager.get_field("Block_Id"))
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



function get_blockNodes(blockID)
    maxBlock = maximum(blockID)
    blockNodes = distribution = [collect(1:maxBlock)]
    for i in 1:maxBlock
        blockNodes[i] = find_indices(blockID, i)
    end
    return blockNodes
end


function solver(params, datamanager)
    blockNodes = init(params, datamanager)
    for block in 1:length(blockNodes)
        datamanager.set_filter(blockNodes[i])
    end
end

function run_solver(blockNodes, datamanager)
    inf = check_inf_or_nan(forces)
    if inf
        @error "Forces are infinite. Time integration is unstable"
        exit()
    end
    #a = (internal_forces + external_forces) / density


end

function solver()
    blockNodes, datamanager = init(params, datamanager)
    run_solver(blockNodes, datamanager)
end

end