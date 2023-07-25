module Thermal_BondBased
include("../../../Support/Parameters/parameter_handling.jl")
export thermal_stretch

function thermal_stretch(params, datamanager)

    alpha = get_element(params, "Alpha")
    temperature = datamanager.get_field("Temperature")
    bondgeom = datamanager.get_field("Deformed Bond Geometry")
    nnodes = datamanager.get_nnodes()
    for iID in 1:nnodes
        bondgeom[iID, :] -= alpha * temperature[iID]
    end
end



end