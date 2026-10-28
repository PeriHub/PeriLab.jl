# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Aniso
using LinearAlgebra: mul!, dot

using ......Data_Manager
using ......PeriLabExceptions: @abort
using ......Helpers: sub_in_place!

export compute_model
export damage_name
export init_model
export fields_for_local_synchronization
"""
    damage_name()

Gives the damage name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: The name of the damage.

Example:
```julia
println(damage_name())
"Critical Energy"
```
"""
function damage_name()
    return "Critical Energy Anisotropic"
end

"""
    compute_model(nodes, damage_parameter, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
[WillbergC2019](@cite), [FosterJT2011](@cite)

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64},
                       damage_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    dof::Int64 = Data_Manager.get_dof()
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    block_ids::NodeScalarField{Int64} = Data_Manager.get_field("Block_Id")
    update_list::NodeScalarField{Bool} = Data_Manager.get_field("Update")
    bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    undeformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Bond Geometry")
    undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
    bond_forces::BondVectorState{Float64} = Data_Manager.get_field("Bond Forces")
    deformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Deformed Bond Geometry",
                                                                     "NP1")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")
    bond_displacements::BondVectorState{Float64} = Data_Manager.get_field("Bond Displacements")
    axis_failed::BondVectorState{Bool} = Data_Manager.get_field("Anisotropic Axis Failure")
    quad_horizons::NodeScalarField{Float64} = Data_Manager.get_field("Quad Horizon")
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
    dissipated_energy::NodeScalarField{Float64} = Data_Manager.get_field("Dissipated Energy")
    inverse_nlist::Vector{Dict{Int64,Int64}} = Data_Manager.get_inverse_nlist()
    rotation_tensor::NodeTensorField{Float64} = Data_Manager.get_field("Rotation Tensor")
    aniso_crit_values::Dict{Int64,Any} = Data_Manager.get_aniso_crit_values()
    tension::Bool = get(damage_parameter, "Only Tension", true)
    aniso_params = get(damage_parameter, "Anisotropic Damage", Dict())
    # "Axis": every material axis fails on its own, the bond keeps the share of the
    # intact axes. "Directional": the bond fails as a whole when its energy exceeds
    # the critical value interpolated for its direction (Ricoeur et al. 2024, Eq. 12).
    directional::Bool = get(aniso_params, "Bond Criterion", "Axis") == "Directional"
    # "Axis Force": a failed axis removes its force component (Correspondence.jl); if the
    # failed axis dominates the bond (w_k >= 0.5) the whole bond fails, so that bonds
    # crossing a crack nearly normal separate the crack faces completely
    axis_force::Bool = get(aniso_params, "Bond Criterion", "Axis") == "Axis Force"
    # stretch based criteria (Zhang & Qiao 2018 transferred to anisotropy); the critical
    # values X, Y, Z are critical stretches s_k instead of energies.
    # "Axis Stretch": axis k fails when the bond stretch exceeds s_k, damage and dominance
    # rule as for "Axis Force". "Directional Stretch": the bond fails as a whole when
    # s^2 > s_min^2 + sum_k (s_k^2 - s_min^2) w_k^(n/2), n = "Exponent".
    axis_stretch::Bool = get(aniso_params, "Bond Criterion", "Axis") == "Axis Stretch"
    directional_stretch::Bool = get(aniso_params, "Bond Criterion", "Axis") ==
                                "Directional Stretch"
    exponent::Float64 = get(aniso_params, "Exponent", 1.0)
    # partly failed bonds fail completely beyond this stretch (test option, off by default)
    stretch_limit::Float64 = get(aniso_params, "Residual Stretch Limit", Inf)

    material_bond = zeros(Float64, dof)
    weights = zeros(Float64, dof)
    relative_disp_buf = zeros(Float64, dof)
    bond_force_delta = zeros(Float64, dof)

    sub_in_place!(bond_displacements, deformed_bond, undeformed_bond)

    for iID in nodes
        rotation = view(rotation_tensor, iID, :, :)
        crit_vals = aniso_crit_values[block_ids[iID]]
        quad_horizon = quad_horizons[iID]

        @inbounds for (jID, neighborID) in enumerate(nlist[iID])
            relative_disp_buf .= bond_displacements[iID][jID]
            norm_displacement = dot(relative_disp_buf, relative_disp_buf)
            if norm_displacement == 0 || (tension &&
                deformed_bond_length[iID][jID] - undeformed_bond_length[iID][jID] < 0)
                continue
            end
            if !haskey(inverse_nlist[neighborID], iID)
                continue
            end

            inID = inverse_nlist[neighborID][iID]
            bond_force_delta .= bond_forces[iID][jID] .- bond_forces[neighborID][inID]
            # same bond energy as in the isotropic Critical Energy model
            bond_energy = 0.25 * abs(dot(bond_force_delta, relative_disp_buf))

            # bond direction in material coordinates (rotation tensor maps material -> global)
            mul!(material_bond, transpose(rotation), deformed_bond[iID][jID])
            # squared direction cosines: sign independent and sum to one, so both half
            # bonds get identical weights and the isotropic limit is the isotropic model
            weights .= (material_bond ./ deformed_bond_length[iID][jID]) .^ 2

            if directional_stretch || axis_stretch
                stretch = deformed_bond_length[iID][jID] /
                          undeformed_bond_length[iID][jID] - 1
                if directional_stretch
                    if stretch^2 >
                       directional_critical_stretch2(crit_vals, weights, exponent)
                        dissipated_energy[iID] += bond_energy * volume[iID] *
                                                  volume[neighborID]
                        bond_damage[iID][jID] = 0.0
                        update_list[iID] = true
                    end
                    continue
                end
                failed = axis_failed[iID][jID]
                new_failure = false
                for k in 1:dof
                    if !failed[k] && stretch > crit_vals[k]
                        failed[k] = true
                        new_failure = true
                    end
                end
                new_failure || continue
                remaining = 1.0
                for k in 1:dof
                    failed[k] || continue
                    remaining -= weights[k]
                    weights[k] >= 0.5 && (remaining = 0.0)
                end
                new_damage = clamp(min(bond_damage[iID][jID], remaining), 0.0, 1.0)
                # bond_energy follows from the damaged bond forces: release the lost share of it
                bond_damage[iID][jID] > 0 &&
                    (dissipated_energy[iID] += bond_energy *
                                               (1 - new_damage / bond_damage[iID][jID]) *
                                               volume[iID] * volume[neighborID])
                bond_damage[iID][jID] = new_damage
                update_list[iID] = true
                continue
            end

            if directional
                if bond_energy >
                   directional_critical_value(crit_vals, weights, exponent) *
                   quad_horizon
                    dissipated_energy[iID] += bond_energy * volume[iID] * volume[neighborID]
                    bond_damage[iID][jID] = 0.0
                    update_list[iID] = true
                end
                continue
            end

            failed = axis_failed[iID][jID]
            new_failure = false
            for k in 1:dof
                if !failed[k] && bond_energy > crit_vals[k] * quad_horizon
                    failed[k] = true
                    new_failure = true
                end
            end
            if !new_failure
                if any(failed) && bond_damage[iID][jID] > 0 &&
                   deformed_bond_length[iID][jID] / undeformed_bond_length[iID][jID] - 1 >
                   stretch_limit
                    dissipated_energy[iID] += bond_energy * volume[iID] * volume[neighborID]
                    bond_damage[iID][jID] = 0.0
                    update_list[iID] = true
                end
                continue
            end

            remaining = 1.0
            for k in 1:dof
                failed[k] || continue
                remaining -= weights[k]
                axis_force && weights[k] >= 0.5 && (remaining = 0.0)
            end
            new_damage = clamp(min(bond_damage[iID][jID], remaining), 0.0, 1.0)
            # bond_energy follows from the damaged bond forces: release the lost share of it
            bond_damage[iID][jID] > 0 &&
                (dissipated_energy[iID] += bond_energy *
                                           (1 - new_damage / bond_damage[iID][jID]) *
                                           volume[iID] * volume[neighborID])
            bond_damage[iID][jID] = new_damage
            update_list[iID] = true
        end
    end
end

"""
    directional_critical_value(crit_vals, weights, exponent)

Critical value of a bond with the squared direction cosines `weights` in material
coordinates, interpolated between the axis values as G = (sum_k w_k G_k^-n)^(-1/n).
For n = 1 this is the harmonic interpolation; in 2D it is Eq. (12) of Ricoeur et al. 2024.
"""
function directional_critical_value(crit_vals, weights::AbstractVector{Float64},
                                    exponent::Float64)
    s = 0.0
    for k in eachindex(weights)
        s += weights[k] * crit_vals[k]^(-exponent)
    end
    return s^(-1 / exponent)
end

"""
    directional_critical_stretch2(crit_stretch, weights, exponent)

Squared critical stretch of a bond with the squared direction cosines `weights`:
s_c^2 = s_min^2 + sum_k (s_k^2 - s_min^2) w_k^(n/2). Bonds along the weakest axis get s_min,
bonds along axis k get s_k; a large n concentrates the higher values on bonds close to the axis.
"""
function directional_critical_stretch2(crit_stretch, weights::AbstractVector{Float64},
                                       exponent::Float64)
    s_min2 = minimum(crit_stretch)^2
    s2 = s_min2
    for k in eachindex(weights)
        s2 += (crit_stretch[k]^2 - s_min2) * weights[k]^(exponent / 2)
    end
    return s2
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    download_from_cores = false
    upload_to_cores = true
    Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

"""
    get_quad_horizon(horizon::Float64, dof::Int64)

Get the quadric of the horizon. The bond energy w is the half bond (one node) share, i.e.
sum_j w_ij V_j is the strain energy density and a bond pair stores 2w. A straight crack is
crossed per unit area by the bonds with int_0^delta 2 xi^2 dxi = 2 delta^3/3 (2D, per thickness),
so the critical half bond energy is w_c = 3 G/(4 delta^3 t). The former 2D value 3/(pi delta^3 t)
dissipated 4/pi G.

# Arguments
- `horizon::Float64`: The horizon of the block.
- `dof::Int64`: The degree of freedom.
- `thickness::Float64`: The thickness of the block.
# Returns
- `quad_horizon::Float64`: The quadric of the horizon.
"""
function get_quad_horizon(horizon::Float64, dof::Int64, thickness::Float64)
    if dof == 2
        return Float64(3 / (4 * horizon^3 * thickness))
    end
    return Float64(4 / (pi * horizon^4))
end

function init_model(nodes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)
    dof = Data_Manager.get_dof()
    quad_horizon = Data_Manager.create_constant_node_scalar_field("Quad Horizon", Float64)
    Data_Manager.create_constant_bond_vector_state("Bond Displacements", Float64, dof)
    # irreversible failure flag per material axis and bond
    Data_Manager.create_constant_bond_vector_state("Anisotropic Axis Failure", Bool, dof)
    Data_Manager.create_constant_node_scalar_field("Dissipated Energy", Float64)
    aniso_params = get(damage_parameter, "Anisotropic Damage", Dict())
    if get(aniso_params, "Bond Criterion", "Axis") == "Axis Force" &&
       get(aniso_params, "Force Filter", true)
        # marker for the correspondence model: a failed axis removes the bond force
        # component along that material axis, the other components stay
        Data_Manager.create_constant_node_scalar_field("Anisotropic Force Filter", Bool)
    end
    horizon = Data_Manager.get_field("Horizon")
    thickness::Float64 = get(damage_parameter, "Thickness", 1)
    for iID in nodes
        quad_horizon[iID] = get_quad_horizon(horizon[iID], dof, thickness)
    end

    if haskey(damage_parameter, "Anisotropic Damage")
        rotation::Bool = Data_Manager.get_rotation()
        if !rotation
            @abort "Anisotropic damage requires Angles field"
            return nothing
        end
    end
end
end
