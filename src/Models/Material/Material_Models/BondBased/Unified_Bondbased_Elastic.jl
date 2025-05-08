# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Unified_Bondbased_Elastic

"""
Based on  [GuanJ2024](@cite)
It is simplified. Rigid body motion is not considered )omega_{ij} = 0 as well as surface correction alpha_i=0 in Equation 16.

"""

#
include("../../Material_Basis.jl")
import .Material_Basis: get_symmetry, apply_pointwise_E
using LoopVectorization
using TimerOutputs
export init_model
export fe_support
export material_name
export compute_model

"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
    return false
end

"""
  init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::Union{SubArray,Vector{Int64}},
                    material_parameter::Dict)
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    constant = datamanager.create_constant_bond_field("Unified Bond Based Constant",
                                                      Float64, 2)
    bb_strain = datamanager.create_constant_node_field("Bond Based Strain", Float64,
                                                       "Matrix", dof)

    bond_length = datamanager.get_field("Bond Length")
    horizon = datamanager.get_field("Horizon")
    symmetry::String = get_symmetry(material_parameter)
    nu = material_parameter["Poisson's Ratio"]
    E = material_parameter["Young's Modulus"]

    # a is not explained -> EQ(41)
    # might be dx
    a = 1
    for iID in nodes
        beta = compute_beta(iID, nu)
        for jID in eachindex(nlist[iID])
            n = bond_length[iID][jID] / horizon[iID]
            I1 = n * sqrt(1 - n * n)
            I2 = n * a * sinh(sqrt(1 / (n * n) - 1))

            if symmetry == "plane stress"
                #constant[iID] = 9 / (pi * horizon[iID]^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
                constant[iID][jID][1],
                constant[iID][jID][2] = compute_constant_2D_pstress(iID,
                                                                    nu,
                                                                    horizon,
                                                                    beta,
                                                                    I1,
                                                                    I2)
            elseif symmetry == "plane strain"
                #constant[iID] = 48 / (5 * pi * horizon[iID]^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
                constant[iID][jID][1],
                constant[iID][jID][2] = compute_constant_2D_pstrain(iID,
                                                                    nu,
                                                                    horizon,
                                                                    beta,
                                                                    I1,
                                                                    I2)
            else
                constant[iID][jID][1] = compute_beta_3D_1(iID, nu) * 12 /
                                        (pi * horizon[iID]^4)
                constant[iID][jID][2] = compute_beta_3D_2(iID, nu) * 12 /
                                        (pi * horizon[iID]^4) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
            end
        end
    end

    return datamanager
end

"""
    material_name()

Returns the name of the material model.
"""
function material_name()
    return "Unified Bond-based Elastic"
end

"""
    compute_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)

Calculate the elastic bond force for each node.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_model(datamanager::Module,
                       nodes::Union{SubArray,Vector{Int64}},
                       material_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    constant = datamanager.get_field("Unified Bond Based Constant")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_damage = datamanager.get_bond_damage("NP1")
    bond_force = datamanager.get_field("Bond Forces")
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    bb_strain = datamanager.get_field("Bond Based Strain")
    E = material_parameter["Young's Modulus"]
    symmetry::String = get_symmetry(material_parameter)

    for iID in nodes
        if any(deformed_bond_length[iID] .== 0)
            @error "Length of bond is zero due to its deformation."
            return nothing
        end
        # Calculate the bond force -> EQ21 for 3D
        # for 2D its 42, 43, 44 and 46
        # -> constants are set up in init (21, 43,44,46)

        @views compute_bond_based_strain(bb_strain[iID, :, :],
                                         deformed_bond[iID],
                                         undeformed_bond[iID],
                                         bond_damage[iID],
                                         nlist[iID],
                                         volume)

        iso_strain = compression_strain(bb_strain[iID, :, :])
        if symmetry != "plane stress" && symmetry != "plane stress"
            compute_bb_force_3D!(bond_force[iID],
                                 0.5 * constant[iID],
                                 bond_damage[iID],
                                 deformed_bond_length[iID],
                                 undeformed_bond_length[iID],
                                 deformed_bond[iID],
                                 iso_strain)
        else
            compute_bb_force_2D!(bond_force[iID],
                                 0.5 * constant[iID],
                                 bond_damage[iID],
                                 deformed_bond_length[iID],
                                 undeformed_bond_length[iID],
                                 deformed_bond[iID],
                                 undeformed_bond[iID],
                                 bond_strain[iID],
                                 iso_strain)
        end
    end
    # might be put in constants
    apply_pointwise_E(nodes, E, bond_force)
    return datamanager
end

"""
based on Equation 16 sum(/deformed_bond-undeformed_bond)/undeformed_bond+ omega)*bond_damage*alpha*V)/sum(bond_damage*alpha*V)
in the paper
    # alpha is volumen correction;
    # omega is rigid body rotation
both is excluded yet
"""

function compute_bond_based_strain(bb_strain,
                                   deformed_bond,
                                   undeformed_bond,
                                   bond_damage,
                                   nlist,
                                   volume)
    reference_volume::Float64 = 0

    reference_volume = compute_reference_volume!(reference_volume, volume, nlist,
                                                 bond_damage)
    # not defined in paper 0/0=0
    @inbounds @fastmath for j in axes(deformed_bond, 1)
        @inbounds @fastmath for k in axes(deformed_bond, 2)
            @inbounds @fastmath for l in axes(undeformed_bond, 2)
                if undeformed_bond[j][l] == 0
                    continue
                end

                bb_strain[k,
                          l] += bond_damage[j] *
                                (deformed_bond[j][k] - undeformed_bond[j][k]) *
                                volume[nlist[j]] /
                                (undeformed_bond[j][l] * reference_volume)
            end
        end
    end
end

function compression_strain(bb_strain)
    iso_strain = zero(eltype(bb_strain))
    @inbounds @fastmath for i in axes(bb_strain, 2)
        iso_strain += bb_strain[i, i]
        #dof += 1
    end
    return iso_strain / 3
end
function compute_reference_volume!(reference_volume::Float64, volume, nlist, bond_damage)
    reference_volume = 0.0
    @inbounds @fastmath for i in eachindex(nlist)
        reference_volume += bond_damage[i] * volume[nlist[i]]
    end
    return reference_volume
end

function compute_deviator(bb_deviator, bb_strain, iso_strain)
    @inbounds @fastmath for j in axes(bb_strain, 1)
        @inbounds @fastmath for i in axes(bb_strain, 2)
            bb_deviator[j, i, i] = bb_strain[j, i, i] - iso_strain[j]
        end
    end
end

function get_beta(beta::Vector{Float64}, nu::Vector{Float64})
    @inbounds @fastmath for i in axes(beta, 1)
        beta[i] = 5 * (1 - 2 * nu[i]) / (2 * (1 + nu[i]))
    end
    return beta
end

function get_beta(beta::Float64, nu::Float64)
    return 5 * (1 - 2 * nu) / (2 * (1 + nu))
end

function compute_bb_force_2D!(bond_force,
                              constant,
                              bond_damage,
                              deformed_bond_length,
                              undeformed_bond_length,
                              deformed_bond) end

function compute_bb_force_3D!(bond_force,
                              constant,
                              bond_damage,
                              deformed_bond_length,
                              undeformed_bond_length,
                              deformed_bond,
                              iso_strain)
    @inbounds @fastmath for i in axes(bond_force, 1)
        @inbounds @fastmath for j in axes(bond_force, 2)
            bond_force[i][j] = bond_damage[i] *
                               (constant[i][1] * iso_strain +
                                constant[i][2] *
                                (deformed_bond_length[i] - undeformed_bond_length[i]) /
                                undeformed_bond_length[i]) *
                               deformed_bond[i][j] / deformed_bond_length[i]
        end
    end
end

function compute_bb_force_2D!(bond_force,
                              constant,
                              bond_damage,
                              deformed_bond_length,
                              undeformed_bond_length,
                              undeformed_bond,
                              deformed_bond,
                              bond_strain,
                              iso_strain)
    @inbounds @fastmath for i in axes(bond_force, 1)
        @inbounds @fastmath for j in axes(bond_force, 2)
            bond_force[i][j] = bond_damage[i] *
                               ((constant[i][2] *
                                 bond_strain[j, j] *
                                 undeformed_bond[i][j] *
                                 undeformed_bond[i][j]) / undeformed_bond_length[i] /
                                undeformed_bond_length[i] +
                                constant[i][2] * iso_strain) *
                               deformed_bond[i][j] / deformed_bond_length[i] # Eq(42)
        end
    end
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

# for 2D its 42, 43, 44 and 46
# -> constants are set up in init (21, 43,44,46)

function compute_beta_3D_1(iID, nu::Union{SubArray,Vector{Float64},Vector{Int64}})
    return 3 * (4 * nu[iID] - 1) / (2 * (1 + nu[iID]))
end
function compute_beta_3D_1(iID, nu::Union{Float64,Int64})
    return 3 * (4 * nu - 1) / (2 * (1 + nu))
end

function compute_beta_3D_2(iID, nu::Union{SubArray,Vector{Float64},Vector{Int64}})
    return 5 * (1 - 2 * nu[iID]) / (2 * (1 + nu[iID]))
end
function compute_beta_3D_2(iID, nu::Union{Float64,Int64})
    return 5 * (1 - 2 * nu) / (2 * (1 + nu))
end

function compute_constant_2D_pstress(iID,
                                     nu::Union{SubArray,Vector{Float64},Vector{Int64}},
                                     horizon,
                                     beta,
                                     I1,
                                     I2)
    c2 = 6 / (pi * horizon[iID] * (1 - nu[iID]))
    Ra = 2 * (1 - nu[iID]) / (1 - 2 * nu[iID]) * beta * I1
    Rb = 6 * nu[iID] * (1 - nu[iID]) / (1 - 2 * nu[iID])^2 * beta * I1 +
         2 * (1 - nu[iID]) / (1 - 2 * nu[iID]) *
         (1 - (1 + nu[iID]) / (1 - 2 * nu[iID]) * beta) *
         I2
    return Ra * c2, Rb * c2
end
function compute_constant_2D_pstress(iID, nu::Union{Float64,Int64}, beta, I1, I2)
    c2 = 6 / (pi * horizon[iID] * (1 - nu))
    Ra = 2 * (1 - nu) / (1 - 2 * nu) * beta * I1
    Rb = 6 * nu * (1 - nu) / (1 - 2 * nu)^2 * beta * I1 +
         2 * (1 - nu) / (1 - 2 * nu) * (1 - (1 + nu) / (1 - 2 * nu) * beta) * I2
    return Ra * c2, Rb * c2
end

function compute_constant_2D_pstrain(iID,
                                     nu::Union{SubArray,Vector{Float64},Vector{Int64}},
                                     horizon,
                                     beta,
                                     I1,
                                     I2)
    c2 = 6 / (pi * horizon[iID] * (1 - 2 * nu[iID]) * (1 + nu[iID]))
    Ra = 2 * (1 + nu[iID]) * beta * I1
    Rb = 2 * (1 + nu[iID]) * (1 - beta) * I2
    return Ra * c2, Rb * c2
end
function compute_constant_2D_pstrain(iID, nu::Union{Float64,Int64}, horizon, beta, I1, I2)
    c2 = 6 / (pi * horizon[iID] * (1 - 2 * nu) * (1 + nu))
    Ra = 2 * (1 + nu) * beta * I1
    Rb = 2 * (1 + nu) * (1 - beta) * I2
    return Ra * c2, Rb * c2
end

function compute_beta(iID, nu::Union{SubArray,Vector{Float64},Vector{Int64}})
    return 5 * (1 - nu[iID]) / (2 * (1 + nu[iID]))
end

function compute_beta(iID, nu::Union{Float64,Int64})
    return 5 * (1 - nu) / (2 * (1 + nu))
end

end
