# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Surface_Correction

using .....Data_Manager
using .....PeriLabExceptions: @abort
using .....InputDeck: SurfaceCorrectionParams
export init_surface_correction
export compute_surface_correction

function compute_surface_correction(nodes,
                                    local_synch,
                                    synchronise_field)
    sc = Data_Manager.get_surface_correction()
    sc === nothing && return
    # SurfaceCorrectionParams allows only "Volume Correction"
    if sc.update
        volumen_correction(nodes, local_synch, synchronise_field)
    end
    return compute_surface_volume_correction(nodes)
end

function compute_surface_volume_correction(nodes)
    nlist = Data_Manager.get_nlist()
    bond_force = Data_Manager.get_field("Bond Forces")
    volume_correction = Data_Manager.get_field("Volume Correction")
    for iID in nodes
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            bond_force[iID][jID][:] .*= volume_correction[iID][jID]
        end
    end
end

function init_surface_correction(sc::Union{Nothing,SurfaceCorrectionParams},
                                 local_synch,
                                 synchronise_field)
    Data_Manager.set_surface_correction(sc)
    sc === nothing && return
    return init_volumen_correction(local_synch, synchronise_field)
end

function init_volumen_correction(local_synch, synchronise_field)
    horizon = Data_Manager.get_field("Horizon")
    vol = Data_Manager.create_constant_node_scalar_field("Reference Volume", Float64)
    volume_correction = Data_Manager.create_constant_bond_scalar_state("Volume Correction",
                                                                       Float64)
    domain_volume = Data_Manager.create_constant_node_scalar_field("Domain Volume", Float64)
    nnodes = Data_Manager.get_nnodes()
    dof = Data_Manager.get_dof()
    for iID in 1:nnodes
        if dof == 3
            vol[iID] = 4 * pi * horizon[iID]^3 / 3
        elseif dof == 2
            vol[iID] = pi * horizon[iID]^2 / 4
        end
    end
    return volumen_correction(1:nnodes, local_synch, synchronise_field)
end
function volumen_correction(nodes, local_synch, synchronise_field)
    nlist = Data_Manager.get_nlist()
    vol = Data_Manager.get_field("Reference Volume")
    volume_correction = Data_Manager.get_field("Volume Correction")
    bond_damage = Data_Manager.get_field("Bond Damage", "NP1")
    volume = Data_Manager.get_field("Volume")
    domain_volume = Data_Manager.get_field("Domain Volume")

    for iID in nodes
        domain_volume[iID] = 0
        for (jID, nID) in enumerate(nlist[iID])
            domain_volume[iID] += volume[nID] * bond_damage[iID][jID]
        end
    end
    Data_Manager.set_local_synch("Surface Correction", "Domain Volume", false, true, 1)
    local_synch("Surface Correction", "upload_to_cores", synchronise_field)
    for iID in nodes
        for (jID, nID) in enumerate(nlist[iID])
            volume_correction[iID][jID] = 2 * vol[iID] /
                                          (domain_volume[nID] + domain_volume[iID])
        end
    end
end
end
