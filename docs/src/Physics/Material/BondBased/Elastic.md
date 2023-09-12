module Elastic_bond_based


function calculate_bond_forces(datamanager)
    """
    Calculates the bond forces based on the bond geometry and deformations.

    ``vec(F) = vec(Deformed Bond Geometry) - vec(Bond Geometry)``

    Parameters:
        - datamanager: DataManager object containing the necessary data.
        
    Returns:
        - Updated datamanager object.
    """
    initial_bond_length = datamanager.get_field("Bond Geometry")
    deformed_bond_length = datamanager.get_field("Deformed Bond Geometry")
    t = datamanager.get_field("Bond Forces")
    nnodes = datamanager.get_nnodes()
    nlist = datamanager.get_nlist()
    dof = datamanager.get_dof()
    bondNumber = 0
    for iID in nnodes
        nneighbors = length(nlist[iID])
        for jID in 1:nneighbors
            bondNumber += 1

            # Calculate the stretch of the bond
            stretch = (deformed_bond_length[bondNumber] - initial_bond_length[bondNumber]) / deformed_bond_length[bondNumber]

            # Calculate the bond force
            t[bondNumber] = 0.5 * (1 - bondDamageNP1[bondNumber]) * stretch * constant
        end
    end
    return datamanager
end
end