# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Surface_correction
export init_surface_correction
export compute_surface_correction

function compute_surface_correction(
    datamanager::Module,
    nodes,
    local_synch,
    synchronise_field,
)
    params = datamanager.get_properties(1, "Surface Correction")
    if isnothing(params["Type"])
        return
    end

    if params["Type"] == "Volume Correction"
        if params["Update"]
            volumen_correction(datamanager, nodes, local_synch, synchronise_field)
        end
        return compute_surface_volume_correction(datamanager, nodes)
    end
end

function compute_surface_volume_correction(datamanager::Module, nodes)
    nlist = datamanager.get_nlist()
    bond_force = datamanager.get_field("Bond Forces")
    volume_correction = datamanager.get_field("Volume Correction")
    for iID in nodes
        for jID in eachindex(nlist[iID])
            bond_force[iID][jID][:] .*= volume_correction[iID][jID]
        end
    end
end

function init_surface_correction(
    datamanager::Module,
    params::Dict,
    local_synch,
    synchronise_field,
)
    # check if surface correction exists
    if !haskey(params, "Surface Correction")
        # if not set it to false
        params["Surface Correction"] = Dict("Type" => nothing)
        datamanager.set_properties("Surface Correction", params["Surface Correction"])
        return datamanager
    end
    # check if type exists; if not its an error
    if !haskey(params["Surface Correction"], "Type")
        @error "Surface Correction needs a Type definition"
        return nothing
    end
    # needed for multi-step, because if type is false its not a valid model
    datamanager.set_properties("Surface Correction", params["Surface Correction"])
    if isnothing(params["Surface Correction"]["Type"])
        return datamanager
    end
    if !haskey(params["Surface Correction"], "Update")
        params["Surface Correction"]["Update"] = false
    end
    if params["Surface Correction"]["Type"] == "Volume Correction"
        return init_volumen_correction(datamanager, local_synch, synchronise_field)
    else
        @error "Type $(params["Surface Correction"]["Type"]) not defined for surface correction."
        return nothing
    end
end

function init_volumen_correction(datamanager, local_synch, synchronise_field)
    horizon = datamanager.get_field("Horizon")
    vol = datamanager.create_constant_node_field("Reference Volume", Float64, 1)
    volume_correction =
        datamanager.create_constant_bond_field("Volume Correction", Float64, 1)
    domain_volume = datamanager.create_constant_node_field("Domain Volume", Float64, 1)
    nnodes = datamanager.get_nnodes()
    dof = datamanager.get_dof()
    for iID = 1:nnodes
        if dof == 3
            vol[iID] = 4 * pi * horizon[iID]^3 / 3
        elseif dof == 2
            vol[iID] = pi * horizon[iID]^2 / 4
        end
    end
    return volumen_correction(datamanager, 1:nnodes, local_synch, synchronise_field)
end
function volumen_correction(datamanager, nodes, local_synch, synchronise_field)


    nlist = datamanager.get_nlist()
    vol = datamanager.get_field("Reference Volume")
    volume_correction = datamanager.get_field("Volume Correction")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    volume = datamanager.get_field("Volume")
    domain_volume = datamanager.get_field("Domain Volume")

    for iID in nodes
        domain_volume[iID] = 0
        for (jID, nID) in enumerate(nlist[iID])
            domain_volume[iID] += volume[nID] * bond_damage[iID][jID]
        end
    end
    datamanager.set_local_synch("Surface Correction", "Domain Volume", false, true, 1)
    local_synch(datamanager, "Surface Correction", "upload_to_cores", synchronise_field)
    for iID in nodes
        for (jID, nID) in enumerate(nlist[iID])
            volume_correction[iID][jID] =
                2 * vol[iID] / (domain_volume[nID] + domain_volume[iID])
        end
    end
    return datamanager
end
end
