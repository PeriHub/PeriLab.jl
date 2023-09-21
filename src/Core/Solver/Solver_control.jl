

module Solver
include("../../Support/helpers.jl")
include("../../Physics/Physics_Factory.jl")
include("../../IO/IO.jl")
include("Verlet.jl")
include("../../Support/Parameters/parameter_handling.jl")
include("../BC_manager.jl")
using .IO
using .Physics
using .Boundary_conditions
using TimerOutputs

function init_bondDamage_and_influence_function(A, B, C)
    for id in eachindex(B)
        B[id] = fill(1.0, size(B[id]))
        C[id] = fill(1.0, size(C[id]))
    end
    return fill(1.0, size(A)), B, C
end

function init(params, datamanager)
    @info "Init Solver"
    # tbd in csv for global vars
    nnodes = datamanager.get_nnodes()
    nslaves = datamanager.get_nslaves()
    allBlockNodes = get_blockNodes(datamanager.get_field("Block_Id"), nnodes + nslaves)
    blockNodes = get_blockNodes(datamanager.get_field("Block_Id"), nnodes)
    density = datamanager.create_constant_node_field("Density", Float32, 1)
    horizon = datamanager.create_constant_node_field("Horizon", Float32, 1)
    density = set_density(params, allBlockNodes, density) # includes the neighbors
    horizon = set_horizon(params, allBlockNodes, horizon) # includes the neighbors
    solver_options = get_solver_options(params)

    omega = datamanager.create_constant_node_field("Influence Function", Float32, 1)
    bondDamageN, bondDamageNP1 = datamanager.create_bond_field("Bond Damage", Float32, 1)
    omega[:], bondDamageN, bondDamageNP1 = init_bondDamage_and_influence_function(omega, bondDamageN, bondDamageNP1)

    if solver_options["Material Models"]
        datamanager = Physics.init_material_model_fields(datamanager)
    end
    if solver_options["Damage Models"]
        datamanager = Physics.init_damage_model_fields(datamanager)
    end
    if solver_options["Thermal Models"]
        datamanager = Physics.init_thermal_model_fields(datamanager)
    end
    if solver_options["Additive Models"]
        datamanager = Physics.init_additive_model_fields(datamanager)
    end

    Physics.read_properties(params, datamanager)
    Physics.init_models(datamanager)
    bcs = Boundary_conditions.init_BCs(params, datamanager)
    if get_solver_name(params) == "Verlet"
        solver_options["Initial Time"], solver_options["dt"], solver_options["nsteps"] = init_Verlet(params, datamanager, blockNodes, solver_options["Material Models"], solver_options["Thermal Models"])
    end
    @info "Finished Init Solver"
    return blockNodes, bcs, datamanager, solver_options
end

function get_blockNodes(blockIDs, nnodes)
    blockNodes = Dict{Int64,Vector{Int64}}()
    for i in unique(blockIDs[1:nnodes])
        blockNodes[i] = find_indices(blockIDs[1:nnodes], i)
    end
    return blockNodes
end

function set_density(params, blockNodes, density)
    for block in eachindex(blockNodes)
        density[blockNodes[block]] .= get_density(params, block)
    end
    return density
end

function set_horizon(params, blockNodes, horizon)
    for block in eachindex(blockNodes)
        horizon[blockNodes[block]] .= get_horizon(params, block)
    end
    return horizon
end


function solver(solver_options, blockNodes, bcs, datamanager, outputs, exos, write_results, to, silent)
    #blockNodes, bcs, datamanager, solver_options = init(params, datamanager)
    # here time steps?
    # run solver -> evaluate; test; and synchro?
    return run_Verlet_solver(solver_options, blockNodes, bcs, datamanager, outputs, exos, write_results, to, silent)

end

function synchronise(comm, datamanager, direction)
    synch_fields = datamanager.get_synch_fields()
    overlap_map = datamanager.get_overlap_map()
    for synch_field in keys(synch_fields)
        if direction == "download_from_cores"
            if synch_fields[synch_field][direction]
                vector = datamanager.get_field(synch_field)
                synch_slaves_to_master(comm, overlap_map, vector, synch_fields[synch_field]["dof"])
            end
        end
        if direction == "upload_to_cores"
            if synch_fields[synch_field][direction]
                vector = datamanager.get_field(synch_field)
                synch_master_to_slaves(comm, overlap_map, vector, synch_fields[synch_field]["dof"])
            end
        end
    end
end

function write_results(exos, dt, outputs, datamanager)
    return IO.write_results(exos, dt, outputs, datamanager)
end

end