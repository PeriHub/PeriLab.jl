# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using CSV
using DataFrames

export get_bc_definitions

"""
    get_bc_definitions(params::Dict)

Get the boundary condition definitions

# Arguments
- `params::Dict`: The parameters
# Returns
- `bcs::Dict{String,Any}`: The boundary conditions
"""
function get_bc_definitions(params::Dict)
    bcs = Dict{String,Any}()
    if haskey(params::Dict, "Boundary Conditions") == false
        return bcs
    end
    for entry in keys(params["Boundary Conditions"])
        bcs[entry] = params["Boundary Conditions"][entry]
    end
    return bcs
end
