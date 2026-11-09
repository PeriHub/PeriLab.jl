# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

# Needs the licensed matrix elements (feature "matrix", local dev mode via .env or a
# license server), like the additive tests.

using PeriLab.Solver_Manager.Matrix_Adapter: init_matrix_system, apply_stiffness!,
                                             matrix_element, field_to_vector!,
                                             vector_to_field!

"Runs the PeriLab setup of `deck` up to the solver initialisation (assembles the serial matrix)."
function init_perilab_matrix(folder::String, deck::String)
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_silent(true)
    input, steps = PeriLab.IO.initialize_data(joinpath(folder, deck), folder,
                                              MPI.COMM_WORLD)
    PeriLab.IO.init_orientations(input.sections.blocks)
    PeriLab.Solver_Manager.init(input, steps[1])
    return nothing
end

"Interleaved (`(i−1)·dof + c`) → component-major (`(c−1)·n + i`) index."
interleaved_to_component(n, dof) = [(c - 1) * n + i for i in 1:n for c in 1:dof]

@testset "ut_perilab_correspondence" begin
    folder = joinpath(@__DIR__, "..", "..", "..", "..", "fullscale_tests",
                      "test_linear_static_matrix_based_solver")
    @testset "matches PeriLab matrix ($deck)" for deck in [
                                                      "linear_static_matrix_based_solver_2D.yaml",
                                                      "linear_static_matrix_based_solver_3D.yaml",
                                                      "linear_static_matrix_based_solver_DCBmodel.yaml"
                                                  ]
        init_perilab_matrix(folder, deck)
        K_perilab = PeriLab.Data_Manager.get_stiffness_matrix()
        sys = init_matrix_system("CPU")
        n = PeriLab.Data_Manager.get_nnodes()
        dof = PeriLab.Data_Manager.get_dof()
        p = interleaved_to_component(n, dof)
        K = to_sparse(sys.K)
        @test size(K) == (n * dof, n * dof)
        @test isapprox(K, K_perilab[p, p]; rtol = 1e-12)

        @testset "apply_stiffness! = force of the serial matrix" begin
            u = rand(MersenneTwister(1), n, dof)
            f = zeros(n, dof)
            apply_stiffness!(f, u, sys)
            @test isapprox(vec(f), -(K_perilab * vec(u)); rtol = 1e-12)
        end

        @testset "field ↔ vector" begin
            u = rand(MersenneTwister(2), n, dof)
            x = zeros(n * dof)
            field_to_vector!(x, u, n)
            @test x[(3 - 1) * dof + dof] == u[3, dof]
            v = zeros(n, dof)
            vector_to_field!(v, x, n)
            @test v == u
        end
    end

    @testset "isolated point" begin
        T = matrix_element("PeriLab Correspondence")
        el = Base.invokelatest(T, 2; zero_energy = true)
        # chain 1 - 2, point 3 without bonds
        data = (nlist_ptr = [1, 2, 3, 3], nlist = [2, 1], volume = [1.0, 1.0, 1.0],
                omega = [1.0, 1.0], bond_damage = [1.0, 1.0],
                bond_geometry = [1.0 -1.0; 0.0 0.0],
                inverse_shape_tensor = repeat([1.0 0.0; 0.0 1.0], 1, 1, 3),
                CB = rand(MersenneTwister(3), 2, 2, 2, 2), zStiff = ones(2, 2, 3))
        K = Base.invokelatest() do
            pattern, K, ke = build_system(el, data, 6)
            compute_elements!(ke, pattern, el, data)
            assemble!(K, pattern, ke)
            Matrix(to_sparse(K))
        end
        @test all(iszero, K[5:6, :])
        @test all(iszero, K[:, 5:6])
        @test any(!iszero, K[1:4, 1:4])
        @test all(isfinite, K * ones(6))
    end
end

@testset "update_system! after bond damage" begin
    folder = joinpath(@__DIR__, "..", "..", "..", "..", "fullscale_tests",
                      "test_linear_static_matrix_based_solver")
    init_perilab_matrix(folder, "linear_static_matrix_based_solver_2D.yaml")
    DM = PeriLab.Data_Manager
    sys = init_matrix_system("CPU")
    n, dof = DM.get_nnodes(), DM.get_dof()
    K_before = to_sparse(sys.K)
    bond_damage = DM.get_bond_damage("NP1")
    update = DM.get_field("Update")
    update .= false
    for (i, k) in ((1, 1), (2, 2), (5, 1))
        bond_damage[i][k] = 0.0
        update[i] = true
    end
    PeriLab.Solver_Manager.Matrix_Adapter.update_system!(sys, update)
    PeriLab.Solver_Manager.Correspondence_matrix_based.rebuild_matrix!(collect(1:n))
    p = interleaved_to_component(n, dof)
    @test !isapprox(to_sparse(sys.K), K_before; rtol = 1e-12)
    @test isapprox(to_sparse(sys.K), DM.get_stiffness_matrix()[p, p]; rtol = 1e-12)
    # nothing flagged: matrix unchanged
    values = copy(sys.K.nzval)
    update .= false
    PeriLab.Solver_Manager.Matrix_Adapter.update_system!(sys, update)
    @test sys.K.nzval == values
end
