# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Solver
include("../../Support/helpers.jl")
include("../../Physics/Physics_Factory.jl")
include("../../Physics/Damage/Damage_Factory.jl")
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

function init(params::Dict, datamanager::Module)
    nnodes = datamanager.get_nnodes()
    nslaves = datamanager.get_nslaves()
    allBlockNodes = get_blockNodes(datamanager.get_field("Block_Id"), nnodes + nslaves)
    blockNodes = get_blockNodes(datamanager.get_field("Block_Id"), nnodes)
    density = datamanager.create_constant_node_field("Density", Float64, 1)
    horizon = datamanager.create_constant_node_field("Horizon", Float64, 1)

    datamanager.create_constant_node_field("Update List", Bool, 1, true)
    density = set_density(params, allBlockNodes, density) # includes the neighbors
    horizon = set_horizon(params, allBlockNodes, horizon) # includes the neighbors
    solver_options = get_solver_options(params)

    datamanager.create_constant_bond_field("Influence Function", Float64, 1, 1)
    datamanager.create_bond_field("Bond Damage", Float64, 1, 1)

    Physics.read_properties(params, datamanager, solver_options["Material Models"])
    datamanager = Physics.init_models(params, datamanager, allBlockNodes, solver_options)
    bcs = Boundary_conditions.init_BCs(params, datamanager)

    if get_solver_name(params) == "Verlet"
        solver_options["Initial Time"], solver_options["dt"], solver_options["nsteps"], solver_options["Numerical Damping"] = Verlet.init_solver(params, datamanager, blockNodes, solver_options["Material Models"], solver_options["Thermal Models"])
    end

    if "Active" in datamanager.get_all_field_keys()
        # can be predefined in mesh. Therefore it should be checked if it is there.
        active = datamanager.get_field("Active")
    else
        active = datamanager.create_constant_node_field("Active", Bool, 1, true)
    end
    @info "Finished Init Solver"
    return blockNodes, bcs, datamanager, solver_options
end

function get_blockNodes(block_ids, nnodes)
    blockNodes = Dict{Int64,Vector{Int64}}()
    for i in unique(block_ids[1:nnodes])
        blockNodes[i] = find_indices(block_ids[1:nnodes], i)
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

function solver(solver_options::Dict{String,Any}, blockNodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Int64,Dict{}}, result_files::Vector{Any}, write_results, to, silent::Bool)

    return Verlet.run_solver(solver_options, blockNodes, bcs, datamanager, outputs, result_files, synchronise_field, write_results, to, silent)

end

function synchronise_field(comm, synch_fields::Dict, overlap_map, get_field, synch_field::String, direction::String)

    if !haskey(synch_fields, synch_field)
        @warn "Field $synch_field does not exists"
        return nothing
    end
    if direction == "download_from_cores"
        if synch_fields[synch_field][direction]
            vector = get_field(synch_field)
            return synch_slaves_to_master(comm, overlap_map, vector, synch_fields[synch_field]["dof"])
        end
        return nothing
    end
    if direction == "upload_to_cores"
        if synch_fields[synch_field][direction]
            vector = get_field(synch_field)
            return synch_master_to_slaves(comm, overlap_map, vector, synch_fields[synch_field]["dof"])
        end
        return nothing
    end
    @error "Wrong direction key word $direction in function synchronise_field; it should be download_from_cores or upload_to_cores"
    return nothing
end

function write_results(result_files, dt, outputs, datamanager)
    return IO.write_results(result_files, dt, outputs, datamanager)
end

end