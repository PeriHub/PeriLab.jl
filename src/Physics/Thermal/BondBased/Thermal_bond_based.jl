module Thermal_BondBased
export thermal_stretch

function thermal_stretch(blockId, datamanager)

    alpha = datamanager.get_element(blockId, "Thermal", "Alpha")
    temperature = datamanager.get_field("Temperature")
    bondgeom = datamanager.get_field("Deformed Bond Geometry")
    nnodes = datamanager.get_nnodes()
    for iID in 1:nnodes
        bondgeom[iID, :] -= alpha * temperature[iID]
    end
end

end