

module Solver
include("../../Support/Parameters/parameter_handling.jl")
include("../../Support/helpers.jl")
include("../../Physics/Physics_Factory.jl")
include("../../IO/IO.jl")
include("../BC_manager.jl")
include("Verlet.jl")

import .IO
import .Physics
import .Boundary_conditions
function init(params, datamanager)

    # tbd in csv for global vars
    blockNodes = get_blockNodes(datamanager.get_field("Block_Id"))
    solver_options = get_solver_options(params)
    if solver_options["Material Models"]
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
    if solver_options["Damage Models"]
        damage = datamanager.create_node_field("Damage", Float32, 1)
        bondDamage = datamanager.create_bond_field("Bond Damage", Float32, 1)
    end
    if solver_options["Thermal Models"]
        temperature = datamanager.create_node_field("Temperature", Float32, 1)
        flow = datamanager.create_node_field("Flow", Float32, dof) # -> check dof

    end
    if solver_options["Additive Models"]
        activated = datamanager.create_node_field("Activated", Bool, 1)
        activated[:] .= false
    end

    Physics.read_properties(params, datamanager)

    bcs = Boundary_conditions.init_BCs(params, datamanager)
    if get_solver_name(params) == "Verlet"
        solver_options["Initial Time"], solver_options["dt"], solver_options["nsteps"] = init_Verlet(params, datamanager, solver_options["Material Models"], solver_options["Thermal Models"])
    end

    return blockNodes, bcs, datamanager, solver_options
end

function get_blockNodes(blockIDs)
    blockNodes = Dict{Int64,Vector{Int64}}()
    for i in unique(blockIDs)
        blockNodes[i] = find_indices(blockIDs, i)
    end
    return blockNodes
end

function solver(params, datamanager, outputs, exos)
    blockNodes, bcs, datamanager, solver_options = init(params, datamanager)
    # here time steps?
    # run solver -> evaluate; test; and synchro?
    Verlet.run_Verlet_solver(solver_options, blockNodes, bcs, datamanage, outputs, exos, write_results)
    return exos
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
function write_results(exos, step, dt, outputs, datamanager)
    IO.write_results(exos, step, dt, outputs, datamanager)
end
end