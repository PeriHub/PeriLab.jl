# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation_Gradient

export compute

"""
    compute(datamanager, nodes)

    Compute the bond deformation gradient.

    # Arguments
    - `datamanager`: Datamanager.
    - `nodes`: List of nodes.
    # Returns
    - `datamanager`: Datamanager.
"""
function compute(datamanager, nodes)
    @warn "Bond_Deformation_Gradient not supported yet."
    return datamanager
end


end