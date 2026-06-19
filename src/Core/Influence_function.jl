# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module Influence_Function
using ..Data_Manager
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

export init_influence_function

"""
    init_influence_function(nodes::AbstractVector{Int64}, params::Dict)

Initializes the influence function field based on the user-specified
parameter "Influence Function". Supports:
  - predefined names (fast path), e.g. "1/xi^2"
  - arbitrary string expressions using the variables `xi`, `xiX`, `xiY`, `xiZ`,
    e.g. "1/xi^2", "exp(-xi/3)", "xiX^2 + xiY^2"
"""
function init_influence_function(nodes::AbstractVector{Int64},
                                 params::Dict)
    if !haskey(params, "Influence Function")
        return
    end

    bond_geometry = Data_Manager.get_field("Bond Geometry")  # BondVectorState
    omega = Data_Manager.get_field("Influence Function")     # BondScalarState
    dof = Data_Manager.get_dof()  # or however dof is determined in your codebase

    expr_str = params["Influence Function"]

    # --- 1) Predefined, fast-path implementations ---
    if expr_str == "1/xi^2"
        influence_function_1_div_xi_squared!(omega, nodes, bond_geometry, dof)
        return
    elseif expr_str == "1"
        influence_function_constant!(omega, nodes, bond_geometry, dof)
        return
    end

    # --- 2) Fallback: generic string expression ---
    influence_function_from_string!(omega, nodes, bond_geometry, dof, expr_str)
    return
end

"""
Fast hard-coded path for 1/xi^2.
bond_geometry[iID][jID] = [xiX, xiY, (xiZ), xi] -> last entry (index dof+1) = bond length.
"""
function influence_function_1_div_xi_squared!(omega, nodes, bond_geometry, dof)
    for iID in nodes
        for jID in eachindex(bond_geometry[iID])
            xi = bond_geometry[iID][jID][dof+1]
            omega[iID][jID] = 1.0 / (xi * xi)
        end
    end
end

function influence_function_constant!(omega, nodes, bond_geometry, dof)
    for iID in nodes
        omega[iID][:] .= 1.0
    end
end

"""
Generic path: builds a compiled function from the given expression once,
then evaluates it per bond.
Variables available in the expression: xi, xiX, xiY, xiZ (xiZ only if dof==3).
"""
function influence_function_from_string!(omega, nodes, bond_geometry, dof, expr_str::String)
    f = build_influence_function(expr_str, dof)

    for iID in nodes
        for jID in eachindex(bond_geometry[iID])
            bond = bond_geometry[iID][jID]
            xiX = bond[1]
            xiY = dof >= 2 ? bond[2] : 0.0
            xiZ = dof == 3 ? bond[3] : 0.0
            xi = bond[dof+1]
            omega[iID][jID] = f(xi, xiX, xiY, xiZ)
        end
    end
end

"""
Parses `expr_str` once into a compiled RuntimeGeneratedFunction
with signature f(xi, xiX, xiY, xiZ).
"""
function build_influence_function(expr_str::String, dof::Int64)
    parsed = Meta.parse(expr_str)

    fn_expr = :(function (xi, xiX, xiY, xiZ)
                    $parsed
                end)

    return @RuntimeGeneratedFunction(fn_expr)
end

end # module Influence_function
