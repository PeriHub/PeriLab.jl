# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module PD_Solid_Elastic
include("../material_basis.jl")
include("./Ordinary/Ordinary.jl")
using TimerOutputs
using StaticArrays
using .Ordinary: compute_weighted_volume, compute_dilatation, calculate_symmetry_params, get_bond_forces
export fe_support
export init_material_model
export material_name
export compute_forces
export init_material_model
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
  init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)

    datamanager.create_constant_node_field("Weighted Volume", Float64, 1)
    datamanager.create_constant_node_field("Dilatation", Float64, 1)

    bond_force_deviatoric_part = datamanager.create_constant_bond_field("Bond Forces Deviatoric", Float64, 1)
    bond_force_isotropic_part = datamanager.create_constant_bond_field("Bond Forces Isotropic", Float64, 1)
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
    compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64, to::TimerOutput)
    
Computes the forces.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `material_parameter::Dict`: The material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64, to::TimerOutput)
    # global dof
    # global nlist
    # global volume
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")

    nneighbors = datamanager.get_field("Number of Neighbors")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_damage = datamanager.get_bond_damage("NP1")
    omega = datamanager.get_field("Influence Function")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    bond_force = datamanager.get_field("Bond Forces")

    bond_force_deviatoric_part = datamanager.get_field("Bond Forces Deviatoric")
    bond_force_isotropic_part = datamanager.get_field("Bond Forces Isotropic")
    # isotropic; deviatoric; all
    weighted_volume = datamanager.get_field("Weighted Volume")
    theta = datamanager.get_field("Dilatation")

    # optimizing, because if no damage it has not to be updated

    @timeit to "Weighted Volume" weighted_volume = compute_weighted_volume(nodes, nlist, undeformed_bond_length, bond_damage, omega, volume)
    @timeit to "Dilatation" theta = compute_dilatation(nodes, nneighbors, nlist, undeformed_bond_length, deformed_bond_length, bond_damage, volume, weighted_volume, omega)
    @timeit to "Bond Forces" bond_force_deviatoric_part, bond_force_isotropic_part = elastic(nodes, dof, undeformed_bond_length, deformed_bond_length, bond_damage, theta, weighted_volume, omega, material_parameter, bond_force_deviatoric_part, bond_force_isotropic_part)
    bond_force = get_bond_forces(nodes, bond_force_deviatoric_part + bond_force_isotropic_part, deformed_bond, deformed_bond_length, bond_force)

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
function elastic(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, undeformed_bond_length::SubArray, deformed_bond_length::SubArray, bond_damage::SubArray, theta::Vector{Float64}, weighted_volume::Vector{Float64}, omega::SubArray, material::Dict, bond_force_deviatoric_part::SubArray, bond_force_isotropic_part::SubArray)

    shear_modulus = material["Shear Modulus"]
    bulk_modulus = material["Bulk Modulus"]
    symmetry::String = get_symmetry(material)
    # kappa::Float64 = 0
    # gamma::Float64 = 0
    # alpha::Float64 = 0
    deviatoric_deformation = @MVector zeros(Float64, dof)

    alpha, gamma, kappa = Ordinary.calculate_symmetry_params(symmetry, shear_modulus, bulk_modulus)

    for iID in nodes
        # Calculate alpha and beta
        if weighted_volume[iID] == 0
            continue
        end
        deviatoric_deformation = deformed_bond_length[iID] .- undeformed_bond_length[iID] - (gamma * theta[iID] / 3) .* undeformed_bond_length[iID]
        if alpha isa Float64
            bond_force_deviatoric_part[iID] = bond_damage[iID] .* omega[iID] .* alpha .* deviatoric_deformation ./ weighted_volume[iID]
            bond_force_isotropic_part[iID] = bond_damage[iID] .* omega[iID] .* kappa .* theta[iID] .* undeformed_bond_length[iID] ./ weighted_volume[iID]
        else
            bond_force_deviatoric_part[iID] = bond_damage[iID] .* omega[iID] .* alpha[iID] .* deviatoric_deformation ./ weighted_volume[iID]
            bond_force_isotropic_part[iID] = bond_damage[iID] .* omega[iID] .* kappa[iID] .* theta[iID] .* undeformed_bond_length[iID] ./ weighted_volume[iID]
        end
    end
    deviatoric_deformation = nothing
    return bond_force_deviatoric_part, bond_force_isotropic_part
end

end