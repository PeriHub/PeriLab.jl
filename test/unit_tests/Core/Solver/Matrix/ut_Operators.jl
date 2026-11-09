# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "ut_Operators" begin
    n = 20
    K = bar_matrix(n)
    Kd = Matrix(to_sparse(K))
    p = [1, 7, 20]
    f_set = setdiff(1:n, p)
    bc = create_dof_set(p, n)
    x = collect(range(-1.0, 2.0; length = n))
    y = zeros(n)

    @testset "ConstrainedOperator = K_ff on free, identity on prescribed" begin
        mul!(y, ConstrainedOperator(K, bc), x)
        expected = copy(x)
        expected[f_set] = Kd[f_set, f_set] * x[f_set]
        @test y ≈ expected rtol=1e-14
    end

    @testset "constrained_rhs!" begin
        u = zeros(n)
        u[p] = [0.0, 3.0, 5.0]
        f = collect(1.0:n)
        b = zeros(n)
        constrained_rhs!(b, K, bc, u, f)
        @test all(iszero, b[p])
        @test b[f_set] ≈ f[f_set] - Kd[f_set, p] * u[p] rtol=1e-14
    end

    @testset "no sub-matrix memory: allocation independent of problem size (README A6.1)" begin
        function allocation(n)
            K = bar_matrix(n)
            bc = create_dof_set([1, n], n)
            A = ConstrainedOperator(K, bc)
            x = ones(n)
            y = zeros(n)
            mul!(y, A, x)
            return @allocated mul!(y, A, x)
        end
        small = allocation(100)
        large = allocation(20_000)
        # KernelAbstractions 0.9 CPU launches allocate a few hundred bytes (task setup),
        # a vector of 20 000 entries would be 160 kB
        @test small < 4096
        @test large < 4096
    end
end
