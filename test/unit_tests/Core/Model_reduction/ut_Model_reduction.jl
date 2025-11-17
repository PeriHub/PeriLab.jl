# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using MPI
using TimerOutputs

include("../../../src/Core/Model_reduction/Model_reduction.jl")

function guyan_reduction(K::AbstractMatrix{Float64}, m::Vector{Int64}, s::Vector{Int64})
    return K[m, m] - K[m, s]/K[s, s]*K[s, m]
end
