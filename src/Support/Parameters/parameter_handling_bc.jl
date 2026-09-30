# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using CSV
using DataFrames

export get_bc_definitions

"""
    get_bc_definitions(params::AbstractDict)

Get the boundary condition definitions

# Arguments
- `params::AbstractDict`: The parameters
# Returns
- `bcs::AbstractDict{String,Any}`: The boundary conditions
"""
function get_bc_definitions(params::AbstractDict)
    bcs = Dict{String,Any}()
    if haskey(params::AbstractDict, "Boundary Conditions") == false
        return bcs
    end
    for entry in keys(params["Boundary Conditions"])
        bcs[entry] = params["Boundary Conditions"][entry]
    end
    return bcs
end
