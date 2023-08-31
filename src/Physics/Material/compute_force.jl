include("./Ordinary/Ordinary.jl")
function compute_forces(datamanager, time)
    """
        tbd 
        Omega function
        bond damage
        
    """
    nnodes = datamanager.get_nnodes()
    forces = datamanager.get_field("Forces", "NP1")
    nlist = datamanager.get_nlist()
    num = datamanager.get_field("Number of Neighbors")
    volume = datamanager.get_field("Volume")
    bond_geometry = datamanager.get_field("Bond Geometry")

    bond_stretch = datamanager.get_field("Bond Stretch", "NP1")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    theta = datamanager.get_field("Theta")
    weighted_volume = datamanager.get_field("Weighted Volume")
    omega = datamanager.get_field("Influence Function")

    omega[:] .= 1 # Prototype
    bond_damage[:] .= 1# Prototype
    bond_force[iID][jID] = elastic(nnodes, bond_geometry, bond_stretch, bond_damage, theta, weighted_volume, omega)

    for iID in 1:nnodes
        t = elastic(nnodes, bond_geometry, bond_stretch, bond_damage, theta, weighted_volume, omega)
        for jID in 1:num[iID]
            for k in 1:dof
                forces[iID][:, k] = -bond_force[iID][jID] * bond_stretch[iID][jID, k] / bond_stretch[iID][jID, end] * volume[nlist[iID][jID]]
                forces[nlist[iID][jID]][:, k] = bond_force[iID][jID] * bond_stretch[iID][jID, k] / bond_stretch[iID][jID, end] * volume[iID]
            end
        end
        return datamanager
    end
end