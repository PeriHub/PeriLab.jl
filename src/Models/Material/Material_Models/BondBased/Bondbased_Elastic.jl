# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Elastic
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
    constant = datamanager.create_constant_node_field("Bond Based Constant", Float64, 1)
    horizon = datamanager.get_field("Horizon")
    symmetry::String = get_symmetry(material_parameter)
    for iID in nodes
        if symmetry == "plane stress"
            constant[iID] = 12 / (2 * (1 - 1 / 3)) / (pi * horizon[iID]^3) # from EQ 2.9 +2.9 D=2 in Handbook of PD
        elseif symmetry == "plane strain"
            constant[iID] = 12 / (2 * (1 - 0.25 + 0.25 * 0.25)) / (pi * horizon[iID]^3) # from EQ 2.12 + 2.9 D=2 in Handbook of PD
        else
            constant[iID] = 18 / (3 - 2 / 3) / (pi * horizon[iID]^4) # from EQ 2.12 D=3 in Handbook of PD
        end
    end
    return datamanager
end


"""
    material_name()

Returns the name of the material model.
"""
function material_name()
    return "Bond-based Elastic"
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

function compute_bb_force!(
    bond_force,
    constant,
    bond_damage,
    deformed_bond_length,
    undeformed_bond_length,
    deformed_bond,
)
    @inbounds @fastmath for i ∈ axes(bond_force, 1)
        @inbounds @fastmath for j ∈ axes(bond_force, 2)
            bond_force[i, j] =
                (
                    constant *
                    bond_damage[i] *
                    (deformed_bond_length[i] - undeformed_bond_length[i]) /
                    undeformed_bond_length[i]
                ) * deformed_bond[i, j] / deformed_bond_length[i]
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
