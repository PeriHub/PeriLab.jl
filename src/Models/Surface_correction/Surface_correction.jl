# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Surface_correction
export init_surface_correction
export compute_surface_correction

function compute_surface_correction(datamanager::Module, nodes)
    params = datamanager.get_properties(1, "Surface Correction")
    if !params["Type"]
        return
    end
    if params["Type"] == "Volume Correction"
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

    if !haskey(params, "Surface Correction")
        params["Surface Correction"] = Dict("Type" => false)
        for block_id in datamanager.get_block_list()
            datamanager.set_properties(
                block_id,
                "Surface Correction",
                params["Surface Correction"],
            )
        end
        return datamanager
    end
    if !haskey(params["Surface Correction"], "Type")
        @error "Surface Correction needs a Type definition"
        return nothing
    end
    for block_id in datamanager.get_block_list()
        datamanager.set_properties(
            block_id,
            "Surface Correction",
            params["Surface Correction"],
        )
    end
    if params["Surface Correction"]["Type"] == "Volume Correction"
        return init_volumen_correction(datamanager, params, local_synch, synchronise_field)
    else
        @error "Type $(params["Surface Correction"]["Type"]) not defined for surface correction."
        return nothing
    end
end


function init_volumen_correction(datamanager, params, local_synch, synchronise_field)

    nnodes = datamanager.get_nnodes()
    nlist = datamanager.get_nlist()
    horizon = datamanager.get_field("Horizon")
    dof = datamanager.get_dof()
    volume_correction =
        datamanager.create_constant_bond_field("Volume Correction", Float64, 1)
    volume = datamanager.get_field("Volume")
    domain_volume = datamanager.create_constant_node_field("Domain Volume", Float64, 1)

    for iID = 1:nnodes
        domain_volume[iID] = 0
        for jID in nlist[iID]
            domain_volume[iID] += volume[jID]
        end
    end
    datamanager.set_local_synch("Surface Correction", "Domain Volume", false, true, 1)
    local_synch(datamanager, "Surface Correction", "upload_to_cores", synchronise_field)
    for iID = 1:nnodes
        if dof == 3
            vol = 4 * pi * horizon[iID]^3 / 3
        elseif dof == 2
            vol = pi * horizon[iID]^2 / 4
        end
        for (jID, nID) in enumerate(nlist[iID])
            volume_correction[iID][jID] =
                2 * vol / (domain_volume[nID] + domain_volume[iID])
        end
    end
    return datamanager
end






end
