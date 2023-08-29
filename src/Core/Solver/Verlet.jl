
using LinearAlgebra
include("../../Support/tools.jl")
include("../../Support/Parameters/parameter_handling.jl")
function compute_thermodynamic_crititical_time_step(datamanager, lambda, Cv)
    """
    critical time step for a thermodynamic problem
    Selda Oterkus, Erdogan Madenci, and Abigail G. Agwai.  Fully coupled peridynamic thermomechanics
    """
    criticalTimeStep = 1.0e50
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    nnodes = datamanager.get_nnodes()
    density = datamanager.get_field("Density")
    bondgeometry = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    nneighbors = datamanager.get_field("Number of Neighbors")

    lambda = matrix_style(lambda)
    eigLam = maximum(eigvals(lambda))

    for iID in 1:nnodes
        denominator = get_cs_denominator(nneighbors[iID], volume[nlist[iID]], bondgeometry[iID][:, dof+1])
        t = density[iID] * Cv / (eigLam * denominator)
        criticalTimeStep = test_timestep(t, criticalTimeStep)
    end
    return sqrt(criticalTimeStep)
end
function get_cs_denominator(nneighbors, volume, bondgeometry)
    denominator = 0.0
    for jID in 1:nneighbors
        denominator += volume[jID] / bondgeometry[jID]
    end
    return denominator
end
function compute_mechanical_crititical_time_step(datamanager, bulkModulus)
    #https://www.osti.gov/servlets/purl/1140383
    # based on bond-based approximation
    criticalTimeStep = 1.0e50
    nnodes = datamanager.get_nnodes()
    dof = datamanager.get_dof()
    nneighbors = datamanager.get_field("Number of Neighbors")
    nlist = datamanager.get_nlist()
    density = datamanager.get_field("Density")
    bondgeometry = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    horizon = datamanager.get_field("Horizon")


    for iID in 1:nnodes
        denominator = get_cs_denominator(nneighbors[iID], volume[nlist[iID]], bondgeometry[iID][:, dof+1])
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
function compute_crititical_time_step(datamanager, mechanical, thermo)
    criticalTimeStep = 1.0e50
    blocks = datamanager.get_block_list()
    for iblock in blocks
        if thermo
            lambda = datamanager.get_property(iblock, "Thermal Model", "Lambda")
            Cv = datamanager.get_property(iblock, "Thermal Model", "Specific Heat Capacity")
            if (lambda != Nothing) && (Cv != Nothing)
                t = compute_thermodynamic_crititical_time_step(datamanager, lambda, Cv)
                criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)
            end
        end
        if mechanical
            bulkModulus = datamanager.get_property(iblock, "Material Model", "Bulk Modulus")
            if (bulkModulus != Nothing)
                t = compute_mechanical_crititical_time_step(datamanager, bulkModulus)
                criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)
            end
        end
    end
    return criticalTimeStep
end


function init_Verlet(params, datamanager, mechanical, thermo)
    @info "======================="
    @info "==== Verlet Solver ===="
    @info "======================="

    initial_time = get_initial_time(params)
    final_time = get_final_time(params)
    safety_factor = get_safety_factor(params)
    dt = get_fixed_dt(params)
    @info "Initial time: " * string(initial_time) * " [s]"
    @info "Final time: " * string(final_time) * " [s]"
    if dt == -1
        dt = compute_crititical_time_step(datamanager, mechanical, thermo)
        @info "Minimal time increment: " * string(dt) * " [s]"
    else
        @info "Fixed time increment: " * string(dt) * " [s]"
    end

    nsteps, dt = get_integration_steps(initial_time, final_time, safety_factor * dt)

    @info "Safety Factor: " * string(safety_factor)
    @info "Time increment: " * string(dt) * " [s]"
    @info "Number of steps: " * string(nsteps)
    return initial_time, dt, nsteps
end
function get_integration_steps(initial_time, end_time, dt)
    nsteps = ceil((end_time - initial_time) / dt)
    dt = (end_time - initial_time) / nsteps
    return Int64(nsteps), dt
end


function run_Verlet_solver(solver_options::Dict{String,Real}, blockNodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Any,Any}, exos::Vector{Any})

    dof = datamanager.get_dof()


    forces = datamanager.get_field("Forces", "NP1")
    density = datamanager.get_field("Density")

    uN = datamanager.get_field("Displacements", "N")
    uNP1 = datamanager.get_field("Displacements", "NP1")
    aN = datamanager.get_field("Acceleration", "N")
    aNP1 = datamanager.get_field("Acceleration", "N")
    vN = datamanager.get_field("Velocity", "N")
    vNP1 = datamanager.get_field("Velocity", "N")
    dt = solver_options["dt"]
    datamanager = apply_bc(bcs, datamanager, time)
    time::Float32 = 0.0
    for idt in 1:solver_options["nsteps"]+1
        # one step more, because of init step (time = 0)
        uNP1 = 2 .* uN + vN .* dt + 0.5 * aN .* (dt * dt)
        datamanager = apply_bc(bcs, datamanager, time)
        for block in 1:eachindex(blockNodes)
            datamanager.set_filter(blockNodes[block])
            #forces = Physics.compute_model(datamanager, dt, time)
        end
        check_inf_or_nan(forces, "Forces")
        aNP1 = forces ./ density # element wise
        IO.write_results(datamanager, idt, time, outputs, exos)
        datamanager.switch_NP1_to_N()
        time += dt
    end
end