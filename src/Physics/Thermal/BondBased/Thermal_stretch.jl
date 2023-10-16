# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_BondBased
export thermal_stretch

function thermal_stretch(blockId, datamanager)

    alpha = datamanager.get_property(blockId, "Thermal", "Alpha")
    temperature = datamanager.get_field("Temperature")
    bondgeom = datamanager.get_field("Deformed Bond Geometry")
    nnodes = datamanager.get_nnodes()
    for iID in nnodes
        bondgeom[iID, :] -= alpha * temperature[iID]
    end
end

end