# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Solver
include("../../Support/helpers.jl")
include("../../Physics/Physics_Factory.jl")
include("../../IO/IO.jl")
include("Verlet.jl")
include("../../Support/Parameters/parameter_handling.jl")
include("../BC_manager.jl")
include("../../MPI_communication/MPI_communication.jl")
using .IO
using .Physics
using .Boundary_conditions
using .Verlet
using TimerOutputs

function init_bondDamage_and_influence_function(A, B, C)
    for id in eachindex(A)
        A[id] = fill(1.0, size(A[id]))
        B[id] = fill(1.0, size(B[id]))
        C[id] = fill(1.0, size(C[id]))
    end
    return A, B, C
end

function init(params::Dict, datamanager::Module)
    @info "Init Solver"
    # tbd in csv for global vars
    nnodes = datamanager.get_nnodes()
    nslaves = datamanager.get_nslaves()
    allBlockNodes = get_blockNodes(datamanager.get_field("Block_Id"), nnodes + nslaves)
    blockNodes = get_blockNodes(datamanager.get_field("Block_Id"), nnodes)
    density = datamanager.create_constant_node_field("Density", Float64, 1)
    horizon = datamanager.create_constant_node_field("Horizon", Float64, 1)
    active = datamanager.create_constant_node_field("Active", Bool, 1)
    active .= true
    update_list = datamanager.create_constant_node_field("Update List", Bool, 1)
    update_list .= true
    density = set_density(params, allBlockNodes, density) # includes the neighbors
    horizon = set_horizon(params, allBlockNodes, horizon) # includes the neighbors
    solver_options = get_solver_options(params)

    omega = datamanager.create_constant_bond_field("Influence Function", Float64, 1)
    bondDamageN, bondDamageNP1 = datamanager.create_bond_field("Bond Damage", Float64, 1)
    omega[:], bondDamageN, bondDamageNP1 = init_bondDamage_and_influence_function(omega, bondDamageN, bondDamageNP1)

    if solver_options["Material Models"]
        datamanager = Physics.init_material_model_fields(datamanager)
    end
    if solver_options["Damage Models"]
        datamanager = Physics.init_damage_model_fields(datamanager)
    end
    if solver_options["Thermal Models"]
        datamanager = Physics.init_thermal_model_fields(datamanager)
        heatCapacity = datamanager.create_constant_node_field("Heat Capacity", Float64, 1)
        heatCapacity = set_heatcapacity(params, allBlockNodes, heatCapacity) # includes the neighbors
    end
    if solver_options["Additive Models"]
        datamanager = Physics.init_additive_model_fields(datamanager)
    end

    Physics.read_properties(params, datamanager)
    Physics.init_models(datamanager)
    bcs = Boundary_conditions.init_BCs(params, datamanager)
    if get_solver_name(params) == "Verlet"
        solver_options["Initial Time"], solver_options["dt"], solver_options["nsteps"], solver_options["Numerical Damping"] = Verlet.init_solver(params, datamanager, blockNodes, solver_options["Material Models"], solver_options["Thermal Models"])
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

function set_density(params::Dict, blockNodes::Dict, density::SubArray)
    for block in eachindex(blockNodes)
        density[blockNodes[block]] .= get_density(params, block)
    end
    return density
end

function set_horizon(params::Dict, blockNodes::Dict, horizon::SubArray)
    for block in eachindex(blockNodes)
        horizon[blockNodes[block]] .= get_horizon(params, block)
    end
    return horizon
end



function set_heatcapacity(params::Dict, blockNodes::Dict, heatCapacity::SubArray)
    for block in eachindex(blockNodes)
        heatCapacity[blockNodes[block]] .= get_heatcapacity(params, block)
    end
    return heatCapacity
end
function solver(solver_options::Dict{String,Any}, blockNodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Int64,Dict{}}, result_files::Vector{Any}, write_results, to, silent::Bool)

    return Verlet.run_solver(solver_options, blockNodes, bcs, datamanager, outputs, result_files, synchronise_field, write_results, to, silent)

end

function synchronise_field(comm, synch_fields, overlap_map, get_field, synch_field::String, direction::String)

    if !haskey(synch_fields, synch_field)
        @warn "Field $synch_field does not exists"
        return
    end
    if direction == "download_from_cores"
        if synch_fields[synch_field][direction]
            vector = get_field(synch_field)
            synch_slaves_to_master(comm, overlap_map, vector, synch_fields[synch_field]["dof"])
        end
    end
    if direction == "upload_to_cores"
        if synch_fields[synch_field][direction]
            vector = get_field(synch_field)
            synch_master_to_slaves(comm, overlap_map, vector, synch_fields[synch_field]["dof"])
        end
    end

end

function write_results(result_files, dt, outputs, datamanager)
    return IO.write_results(result_files, dt, outputs, datamanager)
end

end