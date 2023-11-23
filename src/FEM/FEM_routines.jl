# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using LinearAlgebra
using FastGaussQuadrature
function compute_strain(nodes::Union{Vector{Int64},SubArray}, topology, B, u, strain)

    for iID in nodes
        strain[iID, :] = dot(B[iID, :, :], u[topology[iID, :]])

    end
    return strain
end

function get_weights_and_integration_points(dof::Int64, p::Vector{Int64})

    x = zeros(Float64, dof, maximum(p) + 1)
    w = zeros(Float64, dof, maximum(p) + 1)
    for i in 1:dof
        x[i, 1:p[i]+1], w[i, 1:p[i]+1] = gausslegendre(p[i] + 1)
    end
    return x, w
end