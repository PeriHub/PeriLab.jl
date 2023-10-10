# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_stretch
export compute_damage
export damage_name

function damage_name()
    return "Critical Stretch"
end

function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, time::Float32, dt::Float32)

    bondDamageNP1 = datamanager.get_field("Bond Damage", "NP1")
    bondGeom = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    nneighbors = datamanager.get_field("Number of Neighbors")
    cricital_stretch = damage_parameter["Critical Stretch"]
    for iID in nodes
        for jID in nneighbors[iID]
            stretch = (deformed_bond[iID][jID, end] - bondGeom[iID][jID, end]) / bondGeom[iID][jID, end]
            if stretch > cricital_stretch
                bondDamageNP1[iID][jID] = 0
            end
        end
    end
    return datamanager
end

end