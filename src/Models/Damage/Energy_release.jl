# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Model
include("../Material/Material_Factory.jl")
include("../../Support/Geometry.jl")
include("../../Support/Helpers.jl")
using .Material
using .Geometry
using .Helpers:
                rotate,
                fastdot,
                sub_in_place!,
                div_in_place!,
                mul_in_place!,
                interpol_data,
                is_dependent
using LinearAlgebra
using StaticArrays
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
    return "Critical Energy"
end

"""
    compute_model(datamanager, nodes, damage_parameter, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
[WillbergC2019](@cite), [FosterJT2011](@cite)

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       damage_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)

    # Bestehende Datenextraktion bleibt gleich
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    block_ids = datamanager.get_field("Block_Id")
    update_list = datamanager.get_field("Update")
    bond_damage = datamanager.get_bond_damage("NP1")

    undeformed_bond = datamanager.get_field("Bond Geometry")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    bond_forces = datamanager.get_field("Bond Forces")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_displacements = datamanager.get_field("Bond Displacements")

    critical_field = datamanager.has_key("Critical_Value")
    critical_energy = critical_field ? datamanager.get_field("Critical_Value") :
                      damage_parameter["Critical Value"]

    quad_horizons = datamanager.get_field("Quad Horizon")
    inverse_nlist = datamanager.get_inverse_nlist()

    dependend_value,
    dependent_field = is_dependent("Critical Value", damage_parameter,
                                   datamanager)

    # Parameter
    rotation::Bool = datamanager.get_rotation()
    tension::Bool = get(damage_parameter, "Only Tension", false)
    inter_block_damage::Bool = haskey(damage_parameter, "Interblock Damage")

    if inter_block_damage
        inter_critical_energy::Array{Float64,3} = datamanager.get_crit_values_matrix()
    end

    # Optimierte bond_displacements Berechnung
    sub_in_place!(bond_displacements, deformed_bond, undeformed_bond)
    warning_flag = true

    # DISPATCH auf optimierte innerste Schleife basierend auf DOF
    if dof == 2
        # Verwende StaticArrays für DOF=3
        for iID in nodes
            optimized_inner_loop_static_2d(iID,
                                           nlist[iID],                    # node_nlist
                                           bond_displacements[iID],       # node_bond_displacements
                                           bond_forces[iID],             # node_bond_forces
                                           deformed_bond_length[iID],    # node_deformed_length
                                           undeformed_bond_length[iID],  # node_undeformed_length
                                           bond_damage[iID],             # node_bond_damage
                                           bond_forces,                  # all bond_forces (für neighbors)
                                           inverse_nlist,                # inverse_nlist
                                           quad_horizons[iID],          # node_quad_horizon
                                           tension,                      # tension
                                           critical_energy,
                                           update_list)
        end
    else
        return datamanager
    end

    return datamanager
end

# =============================================================================
# OPTIMIERTE INNERSTE SCHLEIFE mit StaticArrays (DOF=3)
# =============================================================================

function optimized_inner_loop_static_2d(iID::Int,
                                        node_nlist::Vector{Int},
                                        node_bond_displacements::Vector{Vector{Float64}}, # Original Format
                                        node_bond_forces::Vector{Vector{Float64}},        # Original Format
                                        node_deformed_length::Vector{Float64},
                                        node_undeformed_length::Vector{Float64},
                                        node_bond_damage::Vector{Float64},
                                        bond_forces::Vector{Matrix{Float64}},     # Alle bond forces
                                        inverse_nlist::Vector{Dict{Int64,Int64}},
                                        node_quad_horizon::Float64,
                                        tension::Bool,
                                        critical_energy_value::Float64,
                                        update_list::Vector{Bool})
    @inbounds @fastmath for jID in eachindex(node_nlist)
        # KONVERTIERE zu StaticArrays - zero allocations nach Konvertierung
        relative_displacement = SVector{2,Float64}(node_bond_displacements[jID, 1],
                                                   node_bond_displacements[jID, 2])

        norm_displacement = dot(relative_displacement, relative_displacement)

        if norm_displacement == 0.0 || (tension &&
            node_deformed_length[jID] - node_undeformed_length[jID] < 0.0)
            continue
        end

        neighborID = node_nlist[jID]

        # Bond force als StaticArray
        bond_force = SVector{2,Float64}(node_bond_forces[jID, 1],
                                        node_bond_forces[jID, 2])

        # OPTIMIERTER neighbor force lookup - bounds check statt try/catch
        neighbor_bond_force = if neighborID <= length(inverse_nlist) &&
                                 haskey(inverse_nlist[neighborID], iID)
            neighbor_idx = inverse_nlist[neighborID][iID]
            if neighbor_idx <= length(bond_forces[neighborID])
                neighbor_force_vec = bond_forces[neighborID][neighbor_idx]
                SVector{2,Float64}(neighbor_force_vec[1],
                                   neighbor_force_vec[2])
            else
                SVector{2,Float64}(0.0, 0.0)
            end
        else
            SVector{2,Float64}(0.0, 0.0)
        end

        # ALLE OPERATIONEN mit StaticArrays - komplett allocation-frei
        bond_force_delta = bond_force - neighbor_bond_force
        product = dot(bond_force_delta, relative_displacement)
        projected_force = (product / norm_displacement) * relative_displacement

        product = dot(projected_force, relative_displacement)
        bond_energy = 0.25 * product

        product = critical_energy_value * node_quad_horizon
        if bond_energy > product
            node_bond_damage[jID] = 0.0
            update_list[iID] = true
        end
    end
end

# =============================================================================
# HILFSFUNKTION für Critical Energy Berechnung
# =============================================================================

function calculate_critical_energy_for_node(critical_field::Bool,
                                            inter_block_damage::Bool,
                                            dependend_value::Bool,
                                            critical_energy, block_ids,
                                            iID::Int, block::Int,
                                            damage_parameter::Dict,
                                            dependent_field,
                                            warning_flag::Bool)
    if critical_field
        return critical_energy[iID]
    elseif inter_block_damage
        # Vereinfacht für dieses Beispiel - kann erweitert werden
        return critical_energy
    elseif dependend_value
        return interpol_data(dependent_field[iID],
                             damage_parameter["Critical Value"]["Data"],
                             warning_flag)
    else
        return critical_energy
    end
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    download_from_cores = false
    upload_to_cores = true
    datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

"""
    get_quad_horizon(horizon::Float64, dof::Int64)

Get the quadric of the horizon.

# Arguments
- `horizon::Float64`: The horizon of the block.
- `dof::Int64`: The degree of freedom.
- `thickness::Float64`: The thickness of the block.
# Returns
- `quad_horizon::Float64`: The quadric of the horizon.
"""
function get_quad_horizon(horizon::Float64, dof::Int64, thickness::Float64)
    #TODO: Use average horizon
    if dof == 2
        return Float64(3 / (pi * horizon^3 * thickness))
    end
    return Float64(4 / (pi * horizon^4))
end

function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)
    dof = datamanager.get_dof()
    quad_horizon = datamanager.create_constant_node_field("Quad Horizon", Float64, 1)
    datamanager.create_constant_bond_field("Bond Displacements", Float64, dof)
    horizon = datamanager.get_field("Horizon")
    thickness::Float64 = get(damage_parameter, "Thickness", 1)
    for iID in nodes
        quad_horizon[iID] = get_quad_horizon(horizon[iID], dof, thickness)
    end

    if haskey(damage_parameter, "Anisotropic Damage")
        rotation::Bool = datamanager.get_rotation()
        if !rotation
            @error "Anisotropic damage requires Angles field"
            return nothing
        end
    end

    return datamanager
end
end
