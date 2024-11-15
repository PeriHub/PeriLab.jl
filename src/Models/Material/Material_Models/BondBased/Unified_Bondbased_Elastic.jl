# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Unified_Bondbased_Elastic

"""
Based on the paper https://doi.org/10.1016/j.apm.2024.01.015
It is simplified. Rigid body motion is not considered )omega_{ij} = 0 as well as surface correction alpha_i=0 in Equation 16.

"""

#
include("../../material_basis.jl")
using .Material_Basis: get_symmetry, apply_pointwise_E
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
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    material_parameter::Dict,
)
    dof = datamanager.get_dof()
    constant =
        datamanager.create_constant_node_field("Unified Bond Based Constant", Float64, 3)
    bb_strain =
        datamanager.create_constant_node_field("Bond Based Strain", Float64, "Matrix", dof)
    horizon = datamanager.get_field("Horizon")
    symmetry::String = get_symmetry(material_parameter)
    nu = material_parameter["Poisson's Modulus"]
    E = material_parameter["Young's Modulus"]
    for iID in nodes
        if symmetry == "plane stress"
            constant[iID] = 12 / (2 * (1 - 1 / 3)) / (pi * horizon[iID]^3) # from EQ 2.9 +2.9 D=2 in Handbook of PD
        elseif symmetry == "plane strain"
            constant[iID] = 12 / (2 * (1 - 0.25 + 0.25 * 0.25)) / (pi * horizon[iID]^3) # from EQ 2.12 + 2.9 D=2 in Handbook of PD
        else
            # beta in equation (29)
            compute_beta(iID, constant[iID, 1], nu)
            constant[iID, 2] = 18 / (3 - 2 / 3) / (pi * horizon[iID]^4) # from EQ 2.12 D=3 in Handbook of PD

        end
    end
    return datamanager
end


function compute_beta(iID, constant, nu::Union{SubArray,Vector{Float64},Vector{Int64}})
    constant = 5 * (1 - 2 * nu[iID]) / (2 * (1 + nu[iID]))
end
function compute_beta(iID, constant, nu::Union{Float64,Int64})
    constant = 5 * (1 - 2 * nu) / (2 * (1 + nu))
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
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    material_parameter::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)

    constant = datamanager.get_field("Bond Based Constant")

    undeformed_bond_length = datamanager.get_field("Bond Length")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_damage = datamanager.get_bond_damage("NP1")
    bond_force = datamanager.get_field("Bond Forces")

    E = material_parameter["Young's Modulus"]

    for iID in nodes
        if any(deformed_bond_length[iID] .== 0)
            @error "Length of bond is zero due to its deformation."
            return nothing
        end
        # Calculate the bond force
        compute_bb_force!(
            bond_force[iID],
            0.5 * constant[iID],
            bond_damage[iID],
            deformed_bond_length[iID],
            undeformed_bond_length[iID],
            deformed_bond[iID],
        )
        #bond_force[iID] =
        #    (
        #        0.5 .* constant[iID] .* bond_damage[iID] .*
        #        (deformed_bond_length[iID] .- undeformed_bond_length[iID]) ./
        #        undeformed_bond_length[iID]
        #    ) .* deformed_bond[iID] ./ deformed_bond_length[iID]
        #
    end
    # checks if E is scalar or a vector. Is needed for point wise definition
    apply_pointwise_E(nodes, E, bond_force)
    return datamanager
end


"""
based on Equation 16
"""

function compute_bond_based_strain(
    bb_strain,
    deformed_bond,
    undeformed_bond,
    bond_damage,
    nlist,
    volume,
)
    reference_volume::Float64 = 0

    reference_volume =
        compute_reference_volume!(reference_volume, volume, nlist, bond_damage)

    @inbounds @fastmath for i ∈ axes(deformed_bond, 1)
        @inbounds @fastmath for j ∈ axes(deformed_bond, 2)
            @inbounds @fastmath for k ∈ axes(undeformed_bond, 2)
                bb_strain[j, k] +=
                    bond_damage[i] *
                    (deformed_bond[i, j] - undeformed_bond[i, j]) *
                    volume[nlist[i]] / (undeformed_bond[i, k] * reference_volume)
            end

        end
    end

    # alpha ist volumen korrektur; omega ist starrkörperrotation
    #summe(/deformed_bond-undeformed_bond)/undeformed_bond+ omega)*bond_damage*alpha*V)/summe(bond_damage*alpha*V)
end


function compression_strain(bb_strain)
    iso_strain = zero(eltype(bb_strain))
    @inbounds @fastmath for i ∈ axes(bb_strain, 1)
        iso_strain += bb_strain[i, i]
    end
    return iso_strain
end
function compute_reference_volume!(reference_volume::Float64, volume, nlist, bond_damage)

    @inbounds @fastmath for i in eachindex(nlist)
        reference_volume += bond_damage[i] * volume[nlist[i]]
    end
    return reference_volume
end


function eq17()

end
function compute_bb_force!(
    bond_force,
    constant,
    bond_damage,
    deformed_bond_length,
    undeformed_bond_length,
    deformed_bond,
)


    #f = beta*c*x**2/bond_norm**2*bond_length/bond_norm  + (1-beta)*c*bond_length/3*bond_length/bond_norm
    #A = 3*(4*nu-1)/((2*(1+nu)))*constant
    #B = 5*(1-nu)/(2*(1+nu))*constant
    #
    #f_bond_norm = A*epsilon_m + B*stretch
    #
    #- Dehnung bestimmen aus
    #
    #epsxyz = deformed_bond/undeformed_bond


    @inbounds @fastmath for i ∈ axes(bond_force, 1)
        @inbounds @fastmath for j ∈ axes(bond_force, 2)
            bond_force[i][j] =
                (
                    constant *
                    bond_damage[i] *
                    (deformed_bond_length[i] - undeformed_bond_length[i]) /
                    undeformed_bond_length[i]
                ) * deformed_bond[i][j] / deformed_bond_length[i]
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


end
