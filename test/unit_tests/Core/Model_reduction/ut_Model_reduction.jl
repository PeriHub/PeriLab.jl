# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using MPI
using TimerOutputs
using SparseArrays
include("../../../../src/Core/Model_reduction/Model_reduction.jl")

using .Model_reduction: stiffness_reduction

@testset "ut_guyan_reduction" begin
    function K_stiff(np, c, nn, L)
        K=zeros(np, np)
        for iID in 1:np
            for jID in (-nn):nn
                if jID != 0 && iID + jID > 0 && iID + jID < np + 1
                    xi = L*abs(jID)
                    K[iID, iID+jID] -= 0.5 * c[iID] / xi
                    K[iID, iID] += 0.5 * c[iID] / xi
                end
            end
        end
        K[np, :] .= 0
        K[:, np] .= 0
        K[np, np] = 1
        return K
    end

    L = 0.12
    nn = 2
    np = 20
    c = ones(np)
    F = zeros(np)
    F[1] = 1.0
    density = ones(np)
    K = sparse(K_stiff(np, c, nn, L))

    u = K \ F

    m = collect(1:10)
    s = collect(11:20)
    m2 = collect(1:13)
    s2 = collect(13:20)
    K_reduced = sparse(stiffness_reduction(K, m, s))
    K_reduced_overlap = sparse(stiffness_reduction(K, m2, s2))
    u_reduced = K_reduced \ F[m]

    #@test isapprox(u[m], u_reduced)
    u1 = rand(20)
    an = collect(11:20)
    u1[an].=0
    F_test = K'*u1
    F_reduced = K_reduced*u1[m]

    F_reduced2 = K_reduced_overlap*u1[m2]
    F_reduced[9:10] = F_reduced2[9:10]
    F_reduced - F_test[m]
    pd_nodes = 1:3
    @test isapprox(F_reduced, F_test[m])
end
