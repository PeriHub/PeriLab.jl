# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Ordinary

export compute_dilatation
export compute_weighted_volume
export get_bond_forces
export calculate_symmetry_params

"""
    compute_weighted_volume(nodes::Union{SubArray,Vector{Int64}},
                            nlist::Vector{Vector{Int64}},
                            undeformed_bond_length::Vector{Vector{Float64}},
                            bond_damage::Vector{Vector{Float64}},
                            omega::Vector{Vector{Float64}},
                            volume::Vector{Float64})

Compute the weighted volume for each node, [SillingSA2007](@cite). Taken from Peridigm -> but adding the bond_damage; this is missing in Peridigm, but should be there.


# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: Vector of node indices or a subarray representing the indices of the nodes.
- `nlist::Vector{Vector{Int64}}`: Vector representing the neighbor list for each node.
- `undeformed_bond_length::Vector{Vector{Float64}}`: Vector representing the undeformed bonds.
- `bond_damage::Vector{Vector{Float64}}`: Vector representing the bond damage.
- `omega::Vector{Vector{Float64}}`: Vector representing the weights for each bond.
- `volume::Vector{Float64}`: Vector representing the volume for each node.

# Returns
- `weighted_volume::Vector{Float64}`: Vector containing the computed weighted volume for each node.
"""
function compute_weighted_volume(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::Vector{Vector{Int64}},
    undeformed_bond_length::Vector{Vector{Float64}},
    bond_damage::Vector{Vector{Float64}},
    omega::Vector{Vector{Float64}},
    volume::Vector{Float64},
)

    if length(nodes) == 0
        return Float64[]
    end
    weighted_volume::Vector{Float64} = zeros(Float64, maximum(nodes))

    for iID in nodes
        # in Peridigm the weighted volume is for some reason independend from damages
        weighted_volume[iID] = sum(
            omega[iID] .* bond_damage[iID] .* undeformed_bond_length[iID] .*
            undeformed_bond_length[iID] .* volume[nlist[iID]],
        )

    end
    # @. weighted_volume[nodes] = sum(omega[nodes] .* bond_damage[nodes] .* undeformed_bond[nodes] .* bond_damage[nodes] .* undeformed_bond[nodes] .* volume[nlist[nodes]], dims=2)
    return weighted_volume
end

"""
    get_bond_forces(nodes::Union{SubArray,Vector{Int64}},
                    bond_force_length::Union{SubArray,Vector{Vector{Float64}}},
                    deformed_bond::SubArray,
                    deformed_bond_length::SubArray,
                    bond_force::SubArray)

Calculate the forces on the bonds in a peridynamic material.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: Vector of node indices or a subarray representing the indices of the nodes.
- `bond_force_length::Union{SubArray,Vector{Vector{Float64}}}`: Vector of vectors or a subarray representing the desired lengths of the bonds for each node.
- `deformed_bond::Vector{Matrix{Float64}}`: Vector representing the deformed bonds.
- `deformed_bond_length::Vector{Vector{Float64}}`: Vector representing the deformed bond lengths.
- `bond_force::Vector{Matrix{Float64}}`: Vector representing the resulting forces on the bonds.

# Returns
- `bond_force::Vector{Matrix{Float64}}`: Vector containing the resulting forces on the bonds.
"""


function get_bond_forces(
    nodes::Union{SubArray,Vector{Int64}},
    bond_force_length::Union{SubArray,Vector{Vector{Float64}}},
    deformed_bond::Vector{Matrix{Float64}},
    deformed_bond_length::Vector{Vector{Float64}},
    bond_force::Vector{Matrix{Float64}},
)
    for iID in nodes
        if any(deformed_bond_length[iID] .== 0)
            @error "Length of bond is zero due to its deformation."
            return nothing
        else
            bond_force[iID] .=
                bond_force_length[iID] .* deformed_bond[iID] ./ deformed_bond_length[iID]
        end
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
function calculate_symmetry_params(
    symmetry::String,
    shear_modulus::Union{Float64,SubArray,Vector{Float64}},
    bulk_modulus::Union{Float64,SubArray,Vector{Float64}},
)
    three_bulk_modulus = 3 .* bulk_modulus
    # from Peridigm damage model. to be checked with literature
    if symmetry == "plane stress"
        return 8 .* shear_modulus,
        4.0 .* shear_modulus ./ (three_bulk_modulus + 4.0 .* shear_modulus),
        4.0 .* bulk_modulus .* shear_modulus ./ (three_bulk_modulus + 4.0 .* shear_modulus)
    elseif symmetry == "plane strain"
        return 8 .* shear_modulus, 2 / 3, (12.0 .* bulk_modulus - 4.0 .* shear_modulus) ./ 9
    else
        return 15 .* shear_modulus, 1, three_bulk_modulus
    end
end

"""
    compute_dilatation(nodes::Union{SubArray,Vector{Int64}}, nneighbors::Vector{Int64}, nlist::Vector{Vector{Int64}}, undeformed_bond::Vector{Vector{Float64}}, deformed_bond::Vector{Vector{Float64}}, deformed_bond_length::Vector{Vector{Float64}}, bond_damage::Vector{Vector{Float64}}, volume::Vector{Float64}, weighted_volume::Vector{Float64}, omega::Vector{Vector{Float64}})

Calculate the dilatation for each node, [SillingSA2007](@cite).

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
- `nneighbors::Vector{Int64}`: Number of neighbors.
- `nlist::Vector{Vector{Int64}}`: Neighbor list.
- `undeformed_bond_length::Vector{Vector{Float64}}`: Bond geometry.
- `deformed_bond_length::Vector{Vector{Float64}}`: Deformed bond geometry.
- `bond_damage::Vector{Vector{Float64}}`: Bond damage.
- `volume::Vector{Float64}`: Volume.
- `weighted_volume::Vector{Float64}`: Weighted volume.
- `omega::Vector{Vector{Float64}}`: Influence function.
# Returns
- `theta::Vector{Float64}`: Dilatation.
"""
function compute_dilatation(
    nodes::Union{SubArray,Vector{Int64}},
    nneighbors::Vector{Int64},
    nlist::Vector{Vector{Int64}},
    undeformed_bond_length::Vector{Vector{Float64}},
    deformed_bond_length::Vector{Vector{Float64}},
    bond_damage::Vector{Vector{Float64}},
    volume::Vector{Float64},
    weighted_volume::Vector{Float64},
    omega::Vector{Vector{Float64}},
)
    # not optimal, because of many zeros, but simpler, because it avoids reorganization. Part of potential optimization
    if length(nodes) == 0
        return Float64[]
    end
    theta::Vector{Float64} = zeros(Float64, maximum(nodes))
    for iID in nodes
        if weighted_volume[iID] == 0
            # @warn "Weighted volume is zero for local point ID: $iID"
            theta[iID] = 0
            continue
        end
        theta[iID] =
            3.0 * sum(
                omega[iID][jID] *
                bond_damage[iID][jID] *
                undeformed_bond_length[iID][jID] *
                (deformed_bond_length[iID][jID] - undeformed_bond_length[iID][jID]) *
                volume[nlist[iID][jID]] / weighted_volume[iID] for jID = 1:nneighbors[iID]
            )
    end
    return theta
end
end