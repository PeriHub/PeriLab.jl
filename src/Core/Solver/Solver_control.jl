

module Solver
include("../../Support/Parameters/parameter_handling.jl")
include("../../Support/helpers.jl")
include("../../Physics/Physics_Factory.jl")

import .IO
import .Physics

function init(params, datamanager)

    exos = init_write_results(output_filenames, datamanager)
    blockNodes = get_blockNodes(datamanager.get_field("Block_Id"))
    mechanical, thermal, additive = get_solver_options(params)
    if mechanical
        dof = datamanager.get_dof()
        force = datamanager.create_node_field("Force", Float32, dof)
        Y = datamanager.create_node_field("Deformed State", Float32, dof)
        u = datamanager.create_node_field("Displacements", Float32, dof)
        bu = datamanager.create_bond_field("Deformed Bond Geometry", Float32, dof + 1)
        a = datamanager.create_constant_node_field("Acceleration", Float32, dof)
        v = datamanager.create_node_field("Velocity", Float32, dof)
        datamanager.set_synch("Force", true, false)
        datamanager.set_synch("Velocity", false, true)
        datamanager.set_synch("Displacements", false, true)
        datamanager.set_synch("Deformed State", false, true)
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

function get_blockNodes(blockIDs)
    blockNodes = Dict{Int64,Vector{Int64}}()
    for i in unique(blockIDs)
        blockNodes[i] = find_indices(blockIDs, i)
    end
    return blockNodes
end

function solver(params, datamanager)
    blockNodes = init(params, datamanager)
    # here time steps?
    # run solver -> evaluate; test; and synchro?
    run_Verlet_solver(blockNodes, datamanager)
end

function run_Verlet_solver(blockNodes, datamanager)

    dof = datamanager.get_dof()
    for block in 1:length(blockNodes)
        datamanager.set_filter(blockNodes[i])
        "evaluate"
    end
    check_inf_or_nan(forces, "Init Forces")

    #init
    #uNP1 = u0 + v0*dt + 0.5*a0*dt*dt

    uNP1 = 2 * uN + v0N * dt + 0.5 * aN * dt * dt

    #a = forces / density


end

function solver()
    blockNodes, datamanager = init(params, datamanager)
    run_solver(blockNodes, datamanager)
    switch_NP1_to_N()
end
function synchronise(comm)
    synch_fields = datamanager.get_synch_fields()
    overlap_map = datamanager.get_overlap_map()
    for synch_field in keys(synch_fields)
        if synch_field["upload_to_cores"]
            set_overlap_information(comm, datamanager.get_field(synch_field), overlap_map)
        end
    end
end
end