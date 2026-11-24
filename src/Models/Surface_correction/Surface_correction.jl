# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Surface_Correction

using MPI
using .....Data_Manager
using .....Helpers: with_dof
using .....InputDeck: SurfaceCorrectionParams
export init_surface_correction
export compute_surface_correction

function compute_surface_correction(nodes,
                                    local_synch,
                                    synchronise_field)
    sc = Data_Manager.get_surface_correction()::Union{Nothing,SurfaceCorrectionParams}
    sc === nothing && return
    # SurfaceCorrectionParams allows only "Volume Correction"
    if sc.update
        volumen_correction(nodes, local_synch, synchronise_field)
    end
    return compute_surface_volume_correction(nodes)
end

function compute_surface_volume_correction(nodes)
    with_dof(D -> _scale_bond_forces!(nodes, D))
end

function _scale_bond_forces!(nodes, ::Val{D}) where {D}
    bond_force = Data_Manager.get_field("Bond Forces")::BondVectorState{Float64,D}
    volume_correction = Data_Manager.get_field("Volume Correction")::BondScalarState{Float64}
    corrected = Data_Manager.get_field("Surface Corrected")::NodeScalarField{Bool}
    @inbounds for iID in nodes
        corrected[iID] || continue
        force = bond_force[iID]
        factor = volume_correction[iID]
        for jID in eachindex(force)
            force[jID] *= factor[jID]
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

# volume of the full horizon; 2D mesh volumes are areas
analytic_volume(horizon, dof) = dof == 3 ? 4 * pi * horizon^3 / 3 : pi * horizon^2

"""
    init_volumen_correction(local_synch, synchronise_field)

The reference volume of a node is the domain volume of the most embedded node of its
block (scaled by the analytical horizon volume, so the horizon may vary). Interior bonds
get the factor 1; the discretization error of the neighbor sum does not rescale the body.
Correspondence blocks are not corrected: the shape tensor accounts for missing neighbors.
"""
function init_volumen_correction(local_synch, synchronise_field)
    horizon = Data_Manager.get_field("Horizon")::NodeScalarField{Float64}
    block_id = Data_Manager.get_field("Block_Id")::NodeScalarField{Int64}
    ref = Data_Manager.create_constant_node_scalar_field("Reference Volume", Float64)
    Data_Manager.create_constant_bond_scalar_state("Volume Correction", Float64;
                                                   default_value = 1.0)
    domain_volume = Data_Manager.create_constant_node_scalar_field("Domain Volume", Float64)
    corrected = Data_Manager.create_constant_node_scalar_field("Surface Corrected", Bool)
    Data_Manager.set_local_synch("Surface Correction", "Domain Volume", false, true, 1)
    nnodes = Data_Manager.get_nnodes()
    dof = Data_Manager.get_dof()
    owned = 1:nnodes

    compute_domain_volume!(domain_volume, owned)
    # ponytail: max over the block; an irregular mesh with an outlier volume gives interior
    # factors below 1, a quantile would be more robust there
    nblocks = maximum(view(block_id, owned); init = 0)
    mpi = Data_Manager.get_mpi_active()
    mpi && (nblocks = MPI.Allreduce(nblocks, MPI.MAX, Data_Manager.get_comm()))
    fill_ratio = zeros(nblocks)
    for iID in owned
        b = block_id[iID]
        fill_ratio[b] = max(fill_ratio[b],
                            domain_volume[iID] / analytic_volume(horizon[iID], dof))
    end
    mpi && (fill_ratio = MPI.Allreduce(fill_ratio, MPI.MAX, Data_Manager.get_comm()))

    skipped = Set{Int64}()
    for iID in owned
        b = block_id[iID]
        ref[iID] = fill_ratio[b] * analytic_volume(horizon[iID], dof)
        material = Data_Manager.get_block_models(b).material
        corrected[iID] = !(material !== nothing && material.correspondence)
        corrected[iID] || push!(skipped, b)
    end
    isempty(skipped) ||
        @info "Surface correction skips the correspondence blocks $(sort!(collect(skipped))): their shape tensor already accounts for missing neighbors"
    return volumen_correction(owned, local_synch, synchronise_field)
end

function compute_domain_volume!(domain_volume, nodes)
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    bond_damage = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    volume = Data_Manager.get_field("Volume")::NodeScalarField{Float64}
    @inbounds for iID in nodes
        v = 0.0
        for (jID, nID) in enumerate(nlist[iID])
            v += volume[nID] * bond_damage[iID][jID]
        end
        domain_volume[iID] = v
    end
end

function volumen_correction(nodes, local_synch, synchronise_field)
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    ref = Data_Manager.get_field("Reference Volume")::NodeScalarField{Float64}
    volume_correction = Data_Manager.get_field("Volume Correction")::BondScalarState{Float64}
    domain_volume = Data_Manager.get_field("Domain Volume")::NodeScalarField{Float64}

    compute_domain_volume!(domain_volume, nodes)
    local_synch("Surface Correction", "upload_to_cores", synchronise_field)
    @inbounds for iID in nodes
        for (jID, nID) in enumerate(nlist[iID])
            # two isolated nodes (all bonds broken): nothing to correct
            sum_volume = domain_volume[nID] + domain_volume[iID]
            volume_correction[iID][jID] = sum_volume > 0 ? 2 * ref[iID] / sum_volume : 1.0
        end
    end
end
end
