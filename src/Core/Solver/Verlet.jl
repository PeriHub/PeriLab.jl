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
function compute_thermodynamic_critical_time_step(nodes::Vector{Int64}, datamanager, lambda, Cv)

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
function compute_mechanical_critical_time_step(nodes::Vector{Int64}, datamanager::Module, bulkModulus::Float32)
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
    return initial_time, dt, nsteps
end
function get_integration_steps(initial_time, end_time, dt)
    if dt <= 0
        @error "Time step $dt [s] is not valid"
    end
    nsteps = ceil((end_time - initial_time) / dt)
    dt = (end_time - initial_time) / nsteps
    return Int64(nsteps), dt
end


function run_solver(solver_options::Dict{String,Any}, blockNodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Int64,Dict{String,Vector{Any}}}, exos::Vector{Any}, synchronise, write_results, to, silent::Bool)
    @info "Run Verlet Solver"
    comm = datamanager.get_comm()
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
    active = datamanager.get_field("Active")
    update_list = datamanager.get_field("Update List")
    dt::Float32 = solver_options["dt"]
    nsteps::Int64 = solver_options["nsteps"]
    start_time::Float32 = solver_options["Initial Time"]
    step_time::Float32 = 0
    for idt in progress_bar(datamanager.get_rank(), nsteps, silent)
        @timeit to "Verlet" begin
            # one step more, because of init step (time = 0)

            vNP1[find_active(active[1:nnodes]), :] = vN[find_active(active[1:nnodes]), :] + 0.5 * dt .* a[find_active(active[1:nnodes]), :]
            # numerical damping vPtr[i] = vPtr[i] * (1 - numericalDamping);
            uNP1[find_active(active[1:nnodes]), :] = uN[find_active(active[1:nnodes]), :] + dt .* vNP1[find_active(active[1:nnodes]), :]

            @timeit to "apply_bc" Boundary_conditions.apply_bc(bcs, datamanager, step_time)
            defCoorNP1[1:nnodes, :] = coor[1:nnodes, :] + uNP1[1:nnodes, :]
            synchronise(comm, datamanager, "upload_to_cores")
            # synch
            for block in eachindex(blockNodes)
                bNodes = blockNodes[block]
                activeBlockNodes = bNodes[find_active(active[bNodes])]
                @timeit to "compute_models" datamanager = Physics.compute_models(datamanager, activeBlockNodes, block, dt, step_time, solver_options, to)
            end
            synchronise(comm, datamanager, "download_from_cores")
            # synch

            check_inf_or_nan(forces_density, "Forces")

            a[find_active(active[1:nnodes]), :] = forces_density[find_active(active[1:nnodes]), :] ./ density[find_active(active[1:nnodes])] # element wise
            forces[find_active(active[1:nnodes]), :] = forces_density[find_active(active[1:nnodes]), :] .* volume[find_active(active[1:nnodes])]

            exos = write_results(exos, start_time + step_time, outputs, datamanager)
            datamanager.switch_NP1_to_N()
            step_time += dt
            if idt < 10 || nsteps - idt < 10 || idt % ceil(nsteps / 10) == 0
                @info "Step: $idt / $nsteps [$step_time s]"
            end
            update_list .= true
        end
    end
    return exos
end

end