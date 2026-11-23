# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "ut_BiCGStab" begin
    n = 60
    data = bar_1d(n; volume = [0.5 + 0.01i for i in 1:n])   # not symmetric
    el = TestSpring{1}()
    pattern, K, ke = build_system(el, data, n)
    compute_elements!(ke, pattern, el, data)
    assemble!(K, pattern, ke)
    p = [1, n]
    bc = create_dof_set(p, n)
    u = zeros(n)
    u[n] = 100.0
    f = zeros(n)
    f[30] = 10.0
    b = zeros(n)
    constrained_rhs!(b, K, bc, u, f)

    ws = BiCGStabWorkspace(b)
    jacobi!(ws, K, bc)
    A = DistributedOperator(ConstrainedOperator(K, bc), nothing)   # serial: no halo
    x, stats = solve_bicgstab!(ws, A, b, SerialSpace(n); rtol = 1e-12)
    apply_prescribed!(x, bc, u)
    reference = reference_solve(Matrix(to_sparse(K)), p, [0.0, 100.0], f)
    @test stats.solved
    @test x ≈ reference rtol=1e-9
    @test x[p] == [0.0, 100.0]

    # second solve with the same workspace allocates nothing (measured inside a function,
    # globals of the test set would allocate themselves)
    # KernelAbstractions 0.9 CPU launches allocate a few hundred bytes each (task setup);
    # per iteration this stays far below one vector of the 2000-DOF bar (16 kB)
    function solve_allocations(n)
        K = bar_matrix(n)
        bc = create_dof_set([1, n], n)
        u = zeros(n)
        u[n] = 100.0
        b = zeros(n)
        constrained_rhs!(b, K, bc, u, zeros(n))
        ws = BiCGStabWorkspace(b)
        jacobi!(ws, K, bc)
        A = DistributedOperator(ConstrainedOperator(K, bc), nothing)
        _, stats = solve_bicgstab!(ws, A, b, SerialSpace(n); rtol = 1e-12)
        return (@allocated solve_bicgstab!(ws, A, b, SerialSpace(n); rtol = 1e-12)) /
               stats.niter
    end
    @test solve_allocations(2000) < 1024
end

@testset "ut_BiCGStab shifted operator and start vector" begin
    n = 40
    K = bar_matrix(n)
    Kd = Matrix(to_sparse(K))
    p = [1, n]
    free = setdiff(1:n, p)
    bc = create_dof_set(p, n)
    shift = [1.0 + 0.1i for i in 1:n]
    b = zeros(n)
    b[free] = sin.(free)
    S = ShiftedOperator(ConstrainedOperator(K, bc), shift)
    y = zeros(n)
    x = collect(range(-1.0, 1.0; length = n))
    mul!(y, S, x)
    @test y[free] ≈ (Kd[free, free] + Diagonal(shift[free])) * x[free] rtol=1e-13
    @test y[p] == x[p]

    ws = BiCGStabWorkspace(b)
    jacobi!(ws, K, bc; shift = shift)
    @test ws.diag_inv[3] ≈ 1 / (Kd[3, 3] + shift[3])
    A = DistributedOperator(S, nothing)
    u, stats = solve_bicgstab!(ws, A, b, SerialSpace(n); rtol = 1e-13)
    reference = zeros(n)
    reference[free] = (Kd[free, free] + Diagonal(shift[free])) \ b[free]
    @test stats.solved
    @test u ≈ reference rtol=1e-10
    # started at the solution: no iteration
    _, stats = solve_bicgstab!(ws, A, b, SerialSpace(n); rtol = 1e-8, x0 = copy(reference))
    @test stats.niter == 0
end
