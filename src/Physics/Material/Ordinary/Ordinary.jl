
module Elastic
include("../material_basis.jl")
export compute_force
export material_name

function material_name()
    return "PD Solid Elastic"
end

function compute_weighted_volume(nnodes, nneighbors, nlist, bond_geometry, bond_damage, omega, volume, weighted_volume)
    weighted_volume = zeros(nnodes)
    for iID in 1:nnodes
        for jID in 1:nneighbors[iID]
            weighted_volume[iID] += omega[iID][jID] * bond_damage[iID][jID] * bond_geometry[iID][jID, end] * bond_geometry[iID][jID, end] * volume[nlist[iID][jID]]
        end

    end
    return weighted_volume
end

function compute_dilatation(nnodes, nneighbors, bond_geometry, deformed_bond, bond_damage, volume, weighted_volume, omega, theta)
    theta = zeros(nnodes)
    for iID in 1:nnodes
        for jID in 1:nneighbors[iID]
            e = deformed_bond[iID][jID, end] - bond_geometry[iID][jID, end]
            theta[iID] += 3.0 * omega[iID][jID] * bond_damage[iID][jID] * bond_geometry[iID][jID, end] * e * volume[nlist[iID][jID]] / weighted_volume[iID]
        end
    end
    return theta
end

function elastic(nnodes, nneighbors, dof, bond_geometry, deformed_bond, bond_damage, theta, weighted_volume, omega, material, bond_force)

    for iID in 1:nnodes
        alpha = 15.0 * material["Shear Modulus"] / weighted_volume[iID]
        beta = 3.0 * material["Bulk Modulus"] / weighted_volume[iID]
        for jID in 1:nneighbors[iID]
            c1 = omega * theta[iID] * (beta - alpha / 3.0)
            t = bond_damage[iID][jID] * omega[iID][jID] * (c1 * bond_geometry[iID][jID, end] + alpha * deformed_bond[iID][jID, end])
            bond_force[iID][jID, :] = t * deformed_bond[iID][jID, 1:dof] / deformed_bond[iID][jID, end]
        end
    end
    return bond_force
end

function compute_force(datamanager, material, time, dt)
    nnodes = datamanager.get_nnodes()
    nlist = datamanager.get_nlist()
    forces = datamanager.get_field("Forces", "NP1")
    nneighbors = datamanager.get_field("Number of Neighbors")
    deformed_bond = datamanager.get_field("Deformed Bond", "NP1")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    omega = datamanager.get_field("Influence Function")
    volume = datamanager.get_field("Volume")
    bond_geometry = datamanager.get_field("Bond Geometry")
    weighted_volume = datamanager.create_constant_node_field("Weighted Volume", 1)
    theta = datamanager.create_constant_node_field("Dilatation", 1)
    bond_force = datamanager.create_constant_bond_field("Bond Forces", 3)
    # optiming, because if no damage it has not to be updated

    weighted_volume = compute_weighted_volume(nnodes, nneighbors, nlist, bond_geometry, bond_damage, omega, volume, weighted_volume)
    theta = compute_dilatation(nnodes, nneighbors, bond_geometry, deformed_bond, bond_damage, volume, weighted_volume, omega, theta)
    bond_force = elastic(nnodes, nneighbors, dof, bond_geometry, deformed_bond, bond_damage, theta, weighted_volume, omega, material, bond_force)
    forces = distribute_forces(nnodes, nneighbors, nlist, bond_force, volume)
    return datamanager
end

end