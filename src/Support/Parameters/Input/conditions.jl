# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct BoundaryConditionParams
    variable::String = req("Variable")
    node_set::String = req("Node Set"; description = "Node set name; several joined with +")
    value::Union{Float64,String} = req("Value"; description = "Number or expression, e.g. \"10*t\"")
    type::Union{Nothing,String} = opt("Type"; default = nothing,
                                      allowed = ["Initial", "Dirichlet", "Neumann"])
    coordinate::Union{Nothing,String} = opt("Coordinate"; default = nothing)
    step_id::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
end

@params struct ComputeClassParams
    compute_class::String = req("Compute Class")
    variable::String = req("Variable")
    calculation_type::Union{Nothing,String} = opt("Calculation Type"; default = nothing)
    block::Union{Nothing,String} = opt("Block"; default = nothing)
    node_set::Union{Nothing,String} = opt("Node Set"; default = nothing)
    equation::Union{Nothing,String} = opt("Equation"; default = nothing)
    x::Union{Nothing,Float64} = opt("X"; default = nothing)
    y::Union{Nothing,Float64} = opt("Y"; default = nothing)
    z::Union{Nothing,Float64} = opt("Z"; default = nothing)
end
