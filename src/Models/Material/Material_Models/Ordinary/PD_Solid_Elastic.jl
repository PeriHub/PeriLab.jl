# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module PD_Solid_Elastic
include("../../Material_Basis.jl")
using .Material_Basis: get_symmetry
include("./Ordinary.jl")
include("../../../../Support/Helpers.jl")
using .Helpers: add_in_place!
using TimerOutputs
using StaticArrays
using .Ordinary:
                 compute_weighted_volume!, compute_dilatation!, calculate_symmetry_params,
                 get_bond_forces
export fe_support
export init_model
export material_name
export compute_model
export init_model
export fields_for_local_synchronization
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
  init_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    datamanager.create_constant_node_field("Weighted Volume", Float64, 1)
    datamanager.create_constant_node_field("Dilatation", Float64, 1)

    bond_force_deviatoric_part = datamanager.create_constant_bond_field("Bond Forces Deviatoric",
                                                                        Float64, 1)
    bond_force_isotropic_part = datamanager.create_constant_bond_field("Bond Forces Isotropic",
                                                                       Float64, 1)
    return datamanager
end

"""
    material_name()

Returns the name of the material model.
"""
function material_name()
    return "PD Solid Elastic"
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

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict, time::Float64, dt::Float64, to::TimerOutput)

Computes the forces.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: The nodes.
- `material_parameter::Dict`: The material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       material_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    # global dof
    # global nlist
    # global volume
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")

    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_damage = datamanager.get_bond_damage("NP1")
    omega = datamanager.get_field("Influence Function")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    bond_force = datamanager.get_field("Bond Forces")
    temp = datamanager.get_field("Temporary Bond Field")

    bond_force_deviatoric_part = datamanager.get_field("Bond Forces Deviatoric")
    bond_force_isotropic_part = datamanager.get_field("Bond Forces Isotropic")
    # isotropic; deviatoric; all
    weighted_volume = datamanager.get_field("Weighted Volume")
    theta = datamanager.get_field("Dilatation")

    # optimizing, because if no damage it has not to be updated
    # TBD update_list should be used here as in shape_tensor.jl
    @timeit to "Weighted Volume" compute_weighted_volume!(weighted_volume,
                                                          nodes,
                                                          nlist,
                                                          undeformed_bond_length,
                                                          bond_damage,
                                                          omega,
                                                          volume)
    @timeit to "Dilatation" compute_dilatation!(nodes,
                                                nlist,
                                                undeformed_bond_length,
                                                deformed_bond_length,
                                                bond_damage,
                                                volume,
                                                weighted_volume,
                                                omega,
                                                theta)
    @timeit to "elastic" elastic!(nodes,
                                  undeformed_bond_length,
                                  deformed_bond_length,
                                  bond_damage,
                                  theta,
                                  weighted_volume,
                                  omega,
                                  material_parameter,
                                  bond_force_deviatoric_part,
                                  bond_force_isotropic_part)
    add_in_place!(temp, bond_force_deviatoric_part, bond_force_isotropic_part)
    @timeit to "get_bond_forces" bond_force=get_bond_forces(nodes, temp, deformed_bond,
                                                            deformed_bond_length,
                                                            bond_force, temp)

    return datamanager
end

"""
    elastic(nodes, dof, undeformed_bond, deformed_bond, bond_damage, theta, weighted_volume, omega, material, bond_force)

Calculate the elastic bond force for each node.

``F = \\omega \\cdot \\theta \\cdot (\\frac{3K}{V} - \\frac{\\frac{15B}{V}}{3} \\cdot \\zeta + \\alpha \\cdot stretch)`` [WillbergC2023](@cite)
for 3D, plane stress and plane strain it is refered to [BobaruF2016](@cite) page 152; Eq. (6.12); after (6.21) and after (6.23)

# Arguments
- nodes: array of node IDs
- dof: number of degrees of freedom
- undeformed_bond: dictionary of bond geometries for each node
- deformed_bond: dictionary of deformed bond geometries for each node
- bond_damage: dictionary of bond damages for each node
- theta: dictionary of theta values for each node
- weighted_volume: dictionary of weighted volumes for each node
- omega: dictionary of omega values for each node
- material: dictionary of material properties
- bond_force: dictionary to store the calculated bond forces for each node

# Returns
- bond_force: dictionary of calculated bond forces for each node
"""
function elastic!(nodes::AbstractVector{Int64},
                  undeformed_bond_length::Vector{Vector{Float64}},
                  deformed_bond_length::Vector{Vector{Float64}},
                  bond_damage::Vector{Vector{Float64}},
                  theta::Vector{Float64},
                  weighted_volume::Vector{Float64},
                  omega::Vector{Vector{Float64}},
                  material::Dict,
                  bond_force_deviatoric_part::Vector{Vector{Float64}},
                  bond_force_isotropic_part::Vector{Vector{Float64}})
    shear_modulus = material["Shear Modulus"]
    bulk_modulus = material["Bulk Modulus"]

    symmetry::String = get_symmetry(material)
    # kappa::Float64 = 0
    # gamma::Float64 = 0
    # alpha::Float64 = 0
    if shear_modulus isa Float64
        alpha, gamma,
        kappa = Ordinary.calculate_symmetry_params(symmetry, shear_modulus,
                                                   bulk_modulus)
    end

    for iID in nodes
        # Calculate alpha and beta
        if weighted_volume[iID] == 0
            continue
        end
        if !(shear_modulus isa Float64)
            alpha, gamma,
            kappa = Ordinary.calculate_symmetry_params(symmetry,
                                                       shear_modulus[iID],
                                                       bulk_modulus[iID])
        end

        deviatoric_deformation::Float64 = 0.0
        @inbounds @simd for jID in eachindex(deformed_bond_length[iID])
            deviatoric_deformation = deformed_bond_length[iID][jID] -
                                     undeformed_bond_length[iID][jID] -
                                     (gamma * theta[iID] / 3) *
                                     undeformed_bond_length[iID][jID]
            bond_force_deviatoric_part[iID][jID] = bond_damage[iID][jID] * omega[iID][jID] *
                                                   alpha *
                                                   deviatoric_deformation /
                                                   weighted_volume[iID]
            bond_force_isotropic_part[iID][jID] = bond_damage[iID][jID] * omega[iID][jID] *
                                                  kappa *
                                                  theta[iID] *
                                                  undeformed_bond_length[iID][jID] /
                                                  weighted_volume[iID]
        end
    end
end

end
