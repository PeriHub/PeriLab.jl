# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Ordinary

include("../../../../Support/Helpers.jl")
using .Helpers: div_in_place!, mul_in_place!
using LinearAlgebra
using LoopVectorization
export compute_dilatation!
export compute_weighted_volume!
export get_bond_forces
export calculate_symmetry_params

"""
    compute_weighted_volume!(nodes::AbstractVector{Int64},
                            nlist::Vector{Vector{Int64}},
                            undeformed_bond_length::Vector{Vector{Float64}},
                            bond_damage::Vector{Vector{Float64}},
                            omega::Vector{Vector{Float64}},
                            volume::Vector{Float64})

Compute the weighted volume for each node, [SillingSA2007](@cite). Taken from Peridigm -> but adding the bond_damage; this is missing in Peridigm, but should be there.


# Arguments
- `nodes::AbstractVector{Int64}`: Vector of node indices or a subarray representing the indices of the nodes.
- `nlist::Vector{Vector{Int64}}`: Vector representing the neighbor list for each node.
- `undeformed_bond_length::Vector{Vector{Float64}}`: Vector representing the undeformed bonds.
- `bond_damage::Vector{Vector{Float64}}`: Vector representing the bond damage.
- `omega::Vector{Vector{Float64}}`: Vector representing the weights for each bond.
- `volume::Vector{Float64}`: Vector representing the volume for each node.

# Returns
- `weighted_volume::Vector{Float64}`: Vector containing the computed weighted volume for each node.
"""
function compute_weighted_volume!(weighted_volume::Vector{Float64},
                                  nodes::AbstractVector{Int64},
                                  nlist::Vector{Vector{Int64}},
                                  undeformed_bond_length::Vector{Vector{Float64}},
                                  bond_damage::Vector{Vector{Float64}},
                                  omega::Vector{Vector{Float64}},
                                  volume::Vector{Float64})
    @inbounds for iID in nodes
        wv = 0.0  # Explizit Float64 statt zero(eltype(weighted_volume))
        @fastmath for (jID, nID) in enumerate(nlist[iID])
            undeformed_length = undeformed_bond_length[iID][jID]
            wv += omega[iID][jID] *
                  bond_damage[iID][jID] *
                  undeformed_length * undeformed_length *  # Vermeidet doppelte Array-Zugriffe
                  volume[nID]
        end
        weighted_volume[iID] = wv
    end
    return nothing
end

"""
    get_bond_forces(nodes::AbstractVector{Int64},
                    bond_force_length::AbstractVector{<:AbstractVector{Float64}},
                    deformed_bond::SubArray,
                    deformed_bond_length::SubArray,
                    bond_force::SubArray)

Calculate the forces on the bonds in a peridynamic material.

# Arguments
- `nodes::AbstractVector{Int64}`: Vector of node indices or a subarray representing the indices of the nodes.
- `bond_force_length::AbstractVector{<:AbstractVector{Float64}}`: Vector of vectors or a subarray representing the desired lengths of the bonds for each node.
- `deformed_bond::Vector{Matrix{Float64}}`: Vector representing the deformed bonds.
- `deformed_bond_length::Vector{Vector{Float64}}`: Vector representing the deformed bond lengths.
- `bond_force::Vector{Matrix{Float64}}`: Vector representing the resulting forces on the bonds.

# Returns
- `bond_force::Vector{Matrix{Float64}}`: Vector containing the resulting forces on the bonds.
"""

function get_bond_forces(nodes::AbstractVector{Int64},
                         bond_force_length::AbstractVector{<:AbstractVector{Float64}},
                         deformed_bond::Vector{Matrix{Float64}},
                         deformed_bond_length::Vector{Vector{Float64}},
                         bond_force::Vector{Matrix{Float64}},
                         temp::Vector{Vector{Float64}})
    div_in_place!(temp, bond_force_length, deformed_bond_length)
    for iID in nodes
        mul_in_place!(bond_force[iID], deformed_bond[iID], temp[iID])
    end
    return bond_force
end

"""
    calculate_symmetry_params(symmetry::String, shear_modulus::Float64, bulk_modulus::Float64)

Calculate the symmetry parameters based on the given material symmetry. These parameters are defined in [BobaruF2016](@cite) Section 6.3, 6.3.1.1 and 6.3.1.2.

# Arguments
- `symmetry::String`: Symmetry of the material.
- `shear_modulus::Float64`: Shear modulus.
- `bulk_modulus::Float64`: Bulk modulus.

# Returns
- `alpha::Float64`: Alpha parameter.
- `gamma::Float64`: Gamma parameter.
- `kappa::Float64`: Kappa parameter.
"""

function calculate_symmetry_params(symmetry::String,
                                   shear_modulus::Float64,
                                   bulk_modulus::Float64)
    three_bulk_modulus = 3 * bulk_modulus
    # from Peridigm damage model. to be checked with literature
    if symmetry == "plane stress"
        return 8 * shear_modulus,
               4.0 * shear_modulus / (three_bulk_modulus + 4.0 * shear_modulus),
               4.0 * bulk_modulus * shear_modulus /
               (three_bulk_modulus + 4.0 * shear_modulus)
    elseif symmetry == "plane strain"
        return 8 * shear_modulus, 2 / 3, (12.0 * bulk_modulus - 4.0 * shear_modulus) / 9
    else
        return 15 * shear_modulus, 1, three_bulk_modulus
    end
end

"""
    compute_dilatation(nodes::AbstractVector{Int64}, nlist::Vector{Vector{Int64}},
                             undeformed_bond_length::Vector{Vector{T}},
                             deformed_bond_length::Vector{Vector{T}},
                             bond_damage::Vector{Vector{T}},
                             volume::Vector{T},
                             weighted_volume::Vector{T},
                             omega::Vector{Vector{T}},
                             theta::Vector{T}) where {T<:AbstractFloat}

Calculate the dilatation for each node, [SillingSA2007](@cite).

# Arguments
- `nodes::AbstractVector{Int64}`: Nodes.
- `nlist::Vector{Vector{Int64}}`: Neighbor list.
- `undeformed_bond_length::Vector{Vector{T}}`: Bond geometry.
- `deformed_bond_length::Vector{Vector{T}}`: Deformed bond geometry.
- `bond_damage::Vector{Vector{T}}`: Bond damage.
- `volume::Vector{T}`: Volume.
- `weighted_volume::Vector{T}`: Weighted volume.
- `omega::Vector{Vector{T}}`: Influence function.
# Returns
- `theta::Vector{Float64}`: Dilatation.
"""

function compute_dilatation!(nodes::AbstractVector{Int64}, nlist::Vector{Vector{Int64}},
                             undeformed_bond_length::Vector{Vector{T}},
                             deformed_bond_length::Vector{Vector{T}},
                             bond_damage::Vector{Vector{T}},
                             volume::Vector{T},
                             weighted_volume::Vector{T},
                             omega::Vector{Vector{T}},
                             theta::Vector{T}) where {T<:AbstractFloat}
    @inbounds for iID in nodes
        if weighted_volume[iID] == zero(T)
            theta[iID] = zero(T)
            continue
        end

        th = zero(T)
        @fastmath @simd for jID in eachindex(nlist[iID])
            th += omega[iID][jID] *
                  bond_damage[iID][jID] *
                  undeformed_bond_length[iID][jID] *
                  (deformed_bond_length[iID][jID] - undeformed_bond_length[iID][jID]) *
                  volume[nlist[iID][jID]]
        end
        theta[iID] = T(3) * th / weighted_volume[iID]
    end
end

end
