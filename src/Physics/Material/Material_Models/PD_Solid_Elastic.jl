# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module PD_Solid_Elastic
include("../material_basis.jl")
include("./Ordinary/Ordinary.jl")
using TimerOutputs
import .Ordinary
export fe_support
export init_material_model
export material_name
export compute_forces

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
    # global dof
    # global nlist
    # global volume

    # dof = datamanager.get_dof()
    # nlist = datamanager.get_nlist()
    # volume = datamanager.get_field("Volume")
    datamanager.create_constant_node_field("Weighted Volume", Float64, 1)
    datamanager.create_constant_node_field("Dilatation", Float64, 1)

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
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    omega = datamanager.get_field("Influence Function")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    bond_force = datamanager.get_field("Bond Forces")
    weighted_volume = datamanager.get_field("Weighted Volume")
    theta = datamanager.get_field("Dilatation")

    # optimizing, because if no damage it has not to be updated

    @timeit to "Weighted Volume" weighted_volume = Ordinary.compute_weighted_volume(nodes, nlist, undeformed_bond, bond_damage, omega, volume)
    @timeit to "Dilatation" theta = Ordinary.compute_dilatation(nodes, nneighbors, nlist, undeformed_bond, deformed_bond, bond_damage, volume, weighted_volume, omega)
    @timeit to "Bond Forces" bond_force = elastic(nodes, dof, undeformed_bond, deformed_bond, bond_damage, theta, weighted_volume, omega, material_parameter, bond_force)


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
function elastic(nodes, dof, undeformed_bond, deformed_bond, bond_damage, theta, weighted_volume, omega, material, bond_force)
    #tbd
    #shear_factor=Vector{Float64}([0,8,15])
    shear_modulus = material["Shear Modulus"]
    bulk_modulus = material["Bulk Modulus"]

    symmetry::String = get_symmmetry(material)
    kappa::Float64 = 0
    gamma::Float64 = 0
    alpha::Float64 = 0
    deviatoric_deformation = zeros(Float64, dof)

    alpha, gamma, kappa = calculate_symmetry_params(symmetry, shear_modulus, bulk_modulus)

    for iID in nodes
        # Calculate alpha and beta
        if weighted_volume[iID] == 0
            continue
        end

        deviatoric_deformation = deformed_bond[iID][:, end] .- undeformed_bond[iID][:, end] - (gamma * theta[iID] / 3) .* undeformed_bond[iID][:, end]
        t = bond_damage[iID][:] .* omega[iID] .* (kappa .* theta[iID] .* undeformed_bond[iID][:, end] .+ alpha .* deviatoric_deformation) ./ weighted_volume[iID]
        if any(deformed_bond[iID][:, end] .== 0)
            @error "Length of bond is zero due to its deformation."
        else
            @inbounds bond_force[iID][:, 1:dof] .= t .* deformed_bond[iID][:, 1:dof] ./ deformed_bond[iID][:, end]
        end
        # Calculate bond force
        #Ordinary.project_bond_forces()
    end
    deviatoric_deformation = nothing

    return bond_force
end

"""
    calculate_symmetry_params(symmetry::String, shear_modulus::Float64, bulk_modulus::Float64)

Calculate the shear modulus, bulk modulus and three bulk modulus for the given symmetry.

# Arguments
- symmetry: symmetry of the material
- shear_modulus: shear modulus
- bulk_modulus: bulk modulus

# Returns
- alpha: alpha
- gamma: gamma
- kappa: kappa
"""
function calculate_symmetry_params(symmetry::String, shear_modulus::Float64, bulk_modulus::Float64)
    three_bulk_modulus = 3 * bulk_modulus
    # from Peridigm damage model. to be checked with literature
    if symmetry == "plane stress"
        return 8 * shear_modulus, 4.0 * shear_modulus / (three_bulk_modulus + 4.0 * shear_modulus), 4.0 * bulk_modulus * shear_modulus / (three_bulk_modulus + 4.0 * shear_modulus)
    elseif symmetry == "plane strain"
        return 8 * shear_modulus, 2 / 3, (12.0 * bulk_modulus - 4.0 * shear_modulus) / 9
    else
        return 15 * shear_modulus, 1, 3 * bulk_modulus
    end
end

end