
using LinearAlgebra
include("../../Support/tools.jl")
function compute_thermodynamic_crititical_time_step(datamanager, lambda, Cv)
    """
    critical time step for a thermodynamic problem
    Selda Oterkus, Erdogan Madenci, and Abigail G. Agwai.  Fully coupled peridynamic thermomechanics
    """
    criticalTimeStep = 1.0e50
    nnodes = datamanager.get_nnodes()
    density = datamanager.get_field("Density")
    bondgeometry = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()

    lambda = matrix_style(lambda)
    eigLam = maximum(eigvals(lambda))

    for iID in 1:nnodes
        denominator = get_cs_denominator(len(nlist[iID]), volume[nlist[iID]], bondgeometry[iID][:, dof+1])
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
function compute_mechanical_crititical_time_step()
    #https://www.osti.gov/servlets/purl/1140383
    # based on bond-based approximation
    criticalTimeStep = 1.0e50
    nnodes = datamanager.get_nnodes()
    density = datamanager.get_field("Density")
    bondgeometry = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    springConstant = 18.0 * bulkModulus / (pi * horizon * horizon * horizon * horizon)
    for iID in 1:nnodes

        springConstant = 18.0 * bulkModulus / (pi * horizon[iID] * horizon[iID] * horizon[iID] * horizon[iID])
        denominator = get_cs_denominator(len(nlist[iID]), volume[nlist[iID]], bondgeometry[iID][:, dof+1])
        t = 18.0 * bulkModulus * density[iID] / (pi * horizon[iID] * horizon[iID] * horizon[iID] * horizon[iID]) / denominator
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
    for iblocks in nblocks
        if thermo
            lambda = 0
            Cv = 0
            t = compute_thermodynamic_crititical_time_step(datamanager, lambda, CV)

            criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)
        end
        if mechanical


            t = compute_mechanical_crititical_time_step()
            criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)


        end
    end
    return criticalTimeStep
end
