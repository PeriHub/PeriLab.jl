
using LinearAlgebra
include("../../Support/tools.jl")

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
function compute_mechanical_crititical_time_step(bulkModulus)
    #https://www.osti.gov/servlets/purl/1140383
    # based on bond-based approximation
    criticalTimeStep = 1.0e50
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    nnodes = datamanager.get_nnodes()
    density = datamanager.get_field("Density")
    bondgeometry = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    horizon = datamanager.get_field("Horizon")

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
    blocks = datamanager.get_block_list()
    for iblock in blocks
        if thermo
            lambda = datamanager.get_property(iblock, "Thermal Models", "Lambda")
            Cv = datamanager.get_property(iblock, "Thermal Models", "Specific Heat Capacity")
            t = compute_thermodynamic_crititical_time_step(datamanager, lambda, Cv)
            criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)
        end
        if mechanical
            bulkModulus = datamanager.get_property(iblock, "Thermal Models", "Bulk Modulus")
            t = compute_mechanical_crititical_time_step(bulkModulus)
            criticalTimeStep = criticalTimeStep = test_timestep(t, criticalTimeStep)
        end
    end
    return criticalTimeStep
end
