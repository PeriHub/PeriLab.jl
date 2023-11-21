# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module PD_Solid_Elastic
include("../material_basis.jl")
include("./Ordinary/Ordinary.jl")
import .Ordinary
export init_material_model
export material_name
export compute_force

# global dof::Int64
# global nlist::Vector{Vector{Int64}}
# global volume::Vector{Float64}

"""
    init_material_model(datamanager::Module)

    Initializes the material model.

    # Arguments
   - `datamanager::Data_manager`: Datamanager.

    # Returns
   - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module)
    # global dof
    # global nlist
    # global volume

    # dof = datamanager.get_dof()
    # nlist = datamanager.get_nlist()
    # volume = datamanager.get_field("Volume")

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
    compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)
    
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
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)
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

    # optimizing, because if no damage it has not to be updated

    weighted_volume = Ordinary.compute_weighted_volume(nodes, nneighbors, nlist, undeformed_bond, bond_damage, omega, volume)
    theta = Ordinary.compute_dilatation(nodes, nneighbors, nlist, undeformed_bond, deformed_bond, bond_damage, volume, weighted_volume, omega)
    bond_force = elastic(nodes, dof, undeformed_bond, deformed_bond, bond_damage, theta, weighted_volume, omega, material_parameter, bond_force)


    return datamanager
end

"""
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

    symmetry::String = get_symmmetry(material)
    kappa::Float64 = 0
    gamma::Float64 = 0
    alpha::Float64 = 0
    deviatoric_deformation = Vector{Float64}
    for iID in nodes
        # Calculate alpha and beta
        if weighted_volume[iID] == 0
            continue
        end
        # from Peridigm damage model. to be checked with literature
        if symmetry == "plane stress"
            alpha = 8 * material["Shear Modulus"]
            gamma = 4.0 * material["Shear Modulus"] / (3.0 * material["Bulk Modulus"] + 4.0 * material["Shear Modulus"])
            kappa = 4.0 * material["Bulk Modulus"] * material["Shear Modulus"] / (3 * material["Bulk Modulus"] + 4.0 * material["Shear Modulus"])
        elseif symmetry == "plane strain"
            alpha = 8 * material["Shear Modulus"]
            gamma = 2 / 3
            kappa = (12.0 * material["Bulk Modulus"] - 4.0 * material["Shear Modulus"]) / 9
        else
            alpha = 15 * material["Shear Modulus"]
            gamma = 1
            kappa = 3 * material["Bulk Modulus"] # -> Eq. (6.12.) 
        end

        deviatoric_deformation = deformed_bond[iID][:, end] .- undeformed_bond[iID][:, end] - (gamma * theta[iID] / 3) .* undeformed_bond[iID][:, end]
        t = bond_damage[iID][:] .* omega[iID] .* (kappa .* theta[iID] .* undeformed_bond[iID][:, end] .+ alpha .* deviatoric_deformation) ./ weighted_volume[iID]
        if deformed_bond[iID][:, end] == 0
            @error "Length of bond is zero due to its deformation."
        end
        # Calculate bond force
        #Ordinary.project_bond_forces()
        bond_force[iID][:, 1:dof] = t .* deformed_bond[iID][:, 1:dof] ./ deformed_bond[iID][:, end]
    end

    return bond_force
end

end