# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


module Verlet
using LinearAlgebra
using TimerOutputs

include("../../Support/helpers.jl")
include("../../Support/tools.jl")
include("../../MPI_communication/MPI_communication.jl")
include("../../Support/Parameters/parameter_handling.jl")
include("../BC_manager.jl")

include("../../Physics/Physics_Factory.jl")
using .Physics
using .Boundary_conditions

export init_solver
export run_solver

"""
[Oterkus2014](@cite)
"""
function compute_thermodynamic_critical_time_step(nodes::Union{SubArray,Vector{Int64}}, datamanager, lambda, Cv)

    criticalTimeStep = 1.0e50
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    density = datamanager.get_field("Density")
    bondgeometry = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    nneighbors = datamanager.get_field("Number of Neighbors")

    lambda = matrix_style(lambda)
    eigLam = maximum(eigvals(lambda))

    for iID in nodes
        denominator = get_cs_denominator(volume[nlist[iID]], bondgeometry[iID][:, dof+1])
        t = density[iID] * Cv / (eigLam * denominator)
        criticalTimeStep = test_timestep(t, criticalTimeStep)
    end
    return sqrt(criticalTimeStep)
end
function get_cs_denominator(volume, bondgeometry)
    return sum(volume ./ bondgeometry)
end
function compute_mechanical_critical_time_step(nodes::Union{SubArray,Vector{Int64}}, datamanager::Module, bulkModulus::Float64)
    #https://www.osti.gov/servlets/purl/1140383
    # based on bond-based approximation
    criticalTimeStep = 1.0e50
    nlist = datamanager.get_nlist()
    density = datamanager.get_field("Density")
    bondgeometry = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    horizon = datamanager.get_field("Horizon")

    for iID in nodes
        denominator = get_cs_denominator(volume[nlist[iID]], bondgeometry[iID][:, end])

        springConstant = 18.0 * bulkModulus / (pi * horizon[iID] * horizon[iID] * horizon[iID] * horizon[iID])

        t = density[iID] / (denominator * springConstant)
        criticalTimeStep = test_timestep(t, criticalTimeStep)
    end
    return sqrt(2 * criticalTimeStep)
end

function test_timestep(t, criticalTimeStep)
    if t < criticalTimeStep
        criticalTimeStep = t
    end
    return criticalTimeStep
end

function compute_crititical_time_step(datamanager, blockNodes::Dict{Int64,Vector{Int64}}, mechanical, thermo)
    criticalTimeStep = 1.0e50
    for iblock in eachindex(blockNodes)
        if thermo
            lambda = datamanager.get_property(iblock, "Thermal Model", "Lambda")
            Cv = datamanager.get_property(iblock, "Thermal Model", "Specific Heat Capacity")
            if (lambda != Nothing) && (Cv != Nothing)
                t = compute_thermodynamic_critical_time_step(blockNodes[iblock], datamanager, lambda, Cv)
                criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)
            end
        end
        if mechanical
            bulkModulus = datamanager.get_property(iblock, "Material Model", "Bulk Modulus")
            if (bulkModulus != Nothing)
                t = compute_mechanical_critical_time_step(blockNodes[iblock], datamanager, bulkModulus)
                criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)
            end
        end
    end
    return criticalTimeStep
end

function init_solver(params, datamanager, blockNodes, mechanical, thermo)
    @info "======================="
    @info "==== Verlet Solver ===="
    @info "======================="

    initial_time = get_initial_time(params)
    final_time = get_final_time(params)
    safety_factor = get_safety_factor(params)
    dt = get_fixed_dt(params)
    @info "Initial time: " * string(initial_time) * " [s]"
    @info "Final time: " * string(final_time) * " [s]"
    if dt == true
        dt = compute_crititical_time_step(datamanager, blockNodes, mechanical, thermo)
        @info "Minimal time increment: " * string(dt) * " [s]"
    else
        @info "Fixed time increment: " * string(dt) * " [s]"
    end

    nsteps, dt = get_integration_steps(initial_time, final_time, safety_factor * dt)
    comm = datamanager.get_comm()
    dt = find_and_set_core_value_min(comm, dt)
    nsteps = find_and_set_core_value_max(comm, nsteps)

    @info "Safety Factor: " * string(safety_factor)
    @info "Time increment: " * string(dt) * " [s]"
    @info "Number of steps: " * string(nsteps)
    @info "Numerical Damping " * string(get_numerical_damping(params))
    return initial_time, dt, nsteps, get_numerical_damping(params)
end
function get_integration_steps(initial_time, end_time, dt)
    if dt <= 0
        @error "Time step $dt [s] is not valid"
    end
    nsteps = ceil((end_time - initial_time) / dt)
    dt = (end_time - initial_time) / nsteps
    return Int64(nsteps), dt
end


function run_solver(solver_options::Dict{String,Any}, blockNodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Int64,Dict{String,Vector{Any}}}, computes, exos::Vector{Any}, csv_files, synchronise_field, write_results, to, silent::Bool)
    @info "Run Verlet Solver"
    dof = datamanager.get_dof()
    nnodes = datamanager.get_nnodes()
    forces = datamanager.get_field("Forces", "NP1")
    volume = datamanager.get_field("Volume")
    forces_density = datamanager.get_field("Force Densities", "NP1")
    density = datamanager.get_field("Density")
    uN = datamanager.get_field("Displacements", "N")
    uNP1 = datamanager.get_field("Displacements", "NP1")
    coor = datamanager.get_field("Coordinates")
    defCoorN = datamanager.get_field("Deformed Coordinates", "N")
    defCoorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
    vN = datamanager.get_field("Velocity", "N")
    vNP1 = datamanager.get_field("Velocity", "NP1")
    a = datamanager.get_field("Acceleration")

    if solver_options["Thermal Models"]
        flowN = datamanager.get_flow("Thermal Flow", "N")
        flowNP1 = datamanager.get_flow("Thermal Flow", "NP1")
        temperatureN = datamanager.get_flow("Temperature", "NP1")
        temperatureNP1 = datamanager.get_flow("Temperature", "NP1")
        heatCapacity = datamanager.get_field("Heat Capacity")
        deltaT = datamanager.create_constant_node_field("Delta Temperature", Float64, 1)
    end
    active = datamanager.get_field("Active")
    update_list = datamanager.get_field("Update List")

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["nsteps"]
    start_time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    numericalDamping::Float64 = solver_options["Numerical Damping"]
    for idt in progress_bar(datamanager.get_rank(), nsteps, silent)
        @timeit to "Verlet" begin
            # one step more, because of init step (time = 0)
            if solver_options["Material Models"]
                vNP1[find_active(active[1:nnodes]), :] = (1 - numericalDamping) .* vN[find_active(active[1:nnodes]), :] + 0.5 * dt .* a[find_active(active[1:nnodes]), :]

                uNP1[find_active(active[1:nnodes]), :] = uN[find_active(active[1:nnodes]), :] + dt .* vNP1[find_active(active[1:nnodes]), :]
            end
            if solver_options["Thermal Models"]
                temperatureNP1[find_active(active[1:nnodes])] = temperatureN[find_active(active[1:nnodes])] + deltaT[find_active(active[1:nnodes])]
            end
            @timeit to "apply_bc" datamanager = Boundary_conditions.apply_bc(bcs, datamanager, step_time)

            defCoorNP1[find_active(active[1:nnodes]), :] = coor[find_active(active[1:nnodes]), :] + uNP1[find_active(active[1:nnodes]), :]

            datamanager.synch_manager(synchronise_field, "upload_to_cores")
            # synch
            for block in eachindex(blockNodes)
                bn = blockNodes[block]
                active_nodes = bn[find_active(active[bn])]
                @timeit to "compute_models" datamanager = Physics.compute_models(datamanager, active_nodes, block, dt, step_time, solver_options, synchronise_field, to)
            end

            datamanager.synch_manager(synchronise_field, "download_from_cores")
            # synch

            check_inf_or_nan(forces_density, "Forces")
            if solver_options["Material Models"]
                a[find_active(active[1:nnodes]), :] = forces_density[find_active(active[1:nnodes]), :] ./ density[find_active(active[1:nnodes])] # element wise
                forces[find_active(active[1:nnodes]), :] = forces_density[find_active(active[1:nnodes]), :] .* volume[find_active(active[1:nnodes])]
            end
            if solver_options["Thermal Models"]
                deltaT[find_active(active[1:nnodes])] = -flowNP1[find_active(active[1:nnodes])] .* dt ./ (density[find_active(active[1:nnodes])] .* heatCapacity[find_active(active[1:nnodes])])
            end
            exos = write_results(exos, csv_files, start_time + step_time, outputs, computes, datamanager)
            datamanager.switch_NP1_to_N()
            update_list .= true
            step_time += dt
            if idt < 10 || nsteps - idt < 10 || idt % ceil(nsteps / 10) == 0
                @info "Step: $idt / $(nsteps+1) [$step_time s]"
            end

        end
    end
    return exos
end

end