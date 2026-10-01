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

    X_kept, w_kept = @test_logs (:warn,) (:info,) CraigBampton.limit_frequency(X, w, 10.0)
    @test w_kept == w
    @test X_kept == X

    X_kept, w_kept = CraigBampton.limit_frequency(X, w, 0.5)
    @test isempty(w_kept)
    @test size(X_kept) == (2, 0)
end
