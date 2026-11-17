# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
using MPI
using TimerOutputs
using SparseArrays

@testset "ut_guyan_reduction" begin
end

@testset "ut_craig_bampton_limit_frequency" begin
    CraigBampton = PeriLab.Solver_Manager.Matrix_Verlet.Model_reduction.CraigBampton
    # eigenvalues of modes at 1, 2 and 3 Hz
    w = (2 * pi .* [1.0, 2.0, 3.0]) .^ 2
    X = [1.0 2.0 3.0; 4.0 5.0 6.0]

    X_kept, w_kept = CraigBampton.limit_frequency(X, w, 2.5)
    @test w_kept == w[1:2]
    @test X_kept == X[:, 1:2]

    # logs are switched off in runtests.jl, so the warning itself is not checked here
    X_kept, w_kept = CraigBampton.limit_frequency(X, w, 10.0)
    @test w_kept == w
    @test X_kept == X

    X_kept, w_kept = CraigBampton.limit_frequency(X, w, 0.5)
    @test isempty(w_kept)
    @test size(X_kept) == (2, 0)
end

@testset "ut_craig_bampton_cascade" begin
    Model_reduction = PeriLab.Solver_Manager.Matrix_Verlet.Model_reduction
    Cascade = Model_reduction.CraigBampton_cascade
    CraigBampton = Model_reduction.CraigBampton
    LinearAlgebra = PeriLab.Solver_Manager.Matrix_Verlet.Model_reduction.CraigBampton.LinearAlgebra

    # 6 x 6 grid, two degrees of freedom per node, numbered component by component
    side = 6
    nnodes = side^2
    dof = 2
    laplace_1d = spdiagm(-1 => -ones(side - 1), 0 => 2 * ones(side), 1 => -ones(side - 1))
    identity_1d = sparse(1.0LinearAlgebra.I, side, side)
    laplace = kron(laplace_1d, identity_1d) + kron(identity_1d, laplace_1d) +
              0.5 * sparse(1.0LinearAlgebra.I, nnodes, nnodes)
    K = kron(sparse([2.0 0.5; 0.5 2.0]), laplace)
    M_diag = 1.0 .+ 0.1 .* (1:(dof * nnodes)) ./ (dof * nnodes)

    retained_nodes = collect(1:side)
    condensed_nodes = collect((side + 1):nnodes)
    r = vcat(retained_nodes, retained_nodes .+ nnodes)
    c = vcat(condensed_nodes, condensed_nodes .+ nnodes)

    subregions = Cascade.split_subregions(c, dof, 3, K)
    @test length(subregions) == 3
    @test sort(vcat(subregions...)) == sort(c)
    # a node is never split between subregions
    for subregion in subregions
        nodes = filter(<=(nnodes), subregion)
        @test sort(subregion) == sort(vcat(nodes, nodes .+ nnodes))
    end
    # wavefront order: the first subregion lies farthest from the retained rows 1 and
    # the last one next to them
    row(node) = (node - 1) ÷ side + 1
    @test minimum(row.(filter(<=(nnodes), subregions[1]))) >
          maximum(row.(filter(<=(nnodes), subregions[end])))

    # static part: exactly the Schur complement, whatever the partition
    K_reduced, M_reduced = Cascade.reduce_matrices(K, M_diag, r, c, 0; dof = dof,
                                                   n_subregions = 3,
                                                   check_symmetry_sample = 0)
    Kcc = Matrix(K[c, c])
    B = Kcc \ Matrix(K[c, r])
    @test Matrix(K_reduced) ≈ Matrix(K[r, r]) - Matrix(K[r, c]) * B
    @test Matrix(M_reduced) ≈
          LinearAlgebra.Diagonal(M_diag[r]) +
          transpose(B) * LinearAlgebra.Diagonal(M_diag[c]) * B

    reduced_eigenvalues(K, M) = sort(real.(LinearAlgebra.eigvals(Matrix(K), Matrix(M))))

    # one subregion spans the same space as the single-region Craig-Bampton
    K_cascade, M_cascade = Cascade.reduce_matrices(K, M_diag, r, c, 4; dof = dof,
                                                   n_subregions = 1,
                                                   check_symmetry_sample = 0)
    K_single, M_single = CraigBampton.reduce_matrices(K, M_diag, r, c, 4;
                                                      check_symmetry_sample = 0)
    @test size(K_cascade) == (length(r) + 4, length(r) + 4)
    @test reduced_eigenvalues(K_cascade, M_cascade) ≈
          reduced_eigenvalues(K_single, M_single)

    # several subregions: n_modes kept in total, and a Ritz approximation from above
    K_cascade, M_cascade = Cascade.reduce_matrices(K, M_diag, r, c, 5; dof = dof,
                                                   n_subregions = 3,
                                                   check_symmetry_sample = 0)
    @test size(K_cascade) == (length(r) + 5, length(r) + 5)
    exact = reduced_eigenvalues(K, LinearAlgebra.Diagonal(M_diag))
    approximate = reduced_eigenvalues(K_cascade, M_cascade)
    @test all(approximate .>= exact[1:length(approximate)] .* (1 - 1e-10))

    # the maximum frequency drops the modes above it
    frequency_limit = 1.2 * sqrt(exact[1]) / (2 * pi)
    K_limited, _ = Cascade.reduce_matrices(K, M_diag, r, c, 5; dof = dof,
                                           n_subregions = 3, check_symmetry_sample = 0,
                                           max_frequency = 1e6)
    @test size(K_limited, 1) == length(r) + 5
    K_limited, _ = Cascade.reduce_matrices(K, M_diag, r, c, 5; dof = dof,
                                           n_subregions = 3, check_symmetry_sample = 0,
                                           max_frequency = 0.0)
    @test size(K_limited, 1) == length(r)
end
