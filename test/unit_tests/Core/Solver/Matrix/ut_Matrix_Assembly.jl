# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "ut_Matrix_Assembly" begin
    @testset "pattern of two small elements" begin
        rows = VectorOfVectors([[1], [1, 2]])
        cols = VectorOfVectors([[1, 3], [2, 3]])
        pattern = create_assembly_pattern(rows, cols, 3)
        @test pattern.ke_ptr == [1, 3, 7]
        @test pattern.rowptr == [1, 4, 6, 6]
        @test pattern.colval == [1, 2, 3, 2, 3]
        K = create_matrix(pattern)
        ke = Float64[1, 2, 10, 20, 30, 40]
        assemble!(K, pattern, ke)
        # K[1,1]=1, K[1,3]=2+30, K[1,2]=10, K[2,2]=20, K[2,3]=40
        @test Matrix(to_sparse(K)) == [1 10 32; 0 20 40; 0 0 0]
    end

    @testset "elements writing into shared rows (correspondence / IGA-PD like)" begin
        rng = Random.MersenneTwister(1)
        n_dofs = 30
        rows = VectorOfVectors{Int64}()
        cols = VectorOfVectors{Int64}()
        for _ in 1:40
            push!(rows, unique(rand(rng, 1:n_dofs, rand(rng, 1:4))))
            push!(cols, unique(rand(rng, 1:n_dofs, rand(rng, 1:8))))
        end
        pattern = create_assembly_pattern(rows, cols, n_dofs)
        ke = rand(rng, pattern.ke_ptr[end] - 1)
        K = create_matrix(pattern)
        assemble!(K, pattern, ke)

        reference = zeros(n_dofs, n_dofs)
        for e in eachindex(rows)
            n_r = length(rows[e])
            for (jc, c) in enumerate(cols[e]), (ir, r) in enumerate(rows[e])
                reference[r, c] += ke[pattern.ke_ptr[e] + (jc - 1) * n_r + ir - 1]
            end
        end
        @test Matrix(to_sparse(K)) ≈ reference rtol=1e-14
        # pattern stays, values are overwritten on reassembly (README A2.5)
        rowptr = copy(K.rowptr)
        assemble!(K, pattern, 2 .* ke)
        @test K.rowptr == rowptr
        @test Matrix(to_sparse(K)) ≈ 2 .* reference rtol=1e-14
    end

    @testset "invalid input" begin
        @test_throws ArgumentError create_assembly_pattern(VectorOfVectors([[4]]),
                                                           VectorOfVectors([[1]]), 3)
        @test_throws ArgumentError create_assembly_pattern(VectorOfVectors([[1]]),
                                                           VectorOfVectors{Int64}(), 3)
    end
end

@testset "status mask (inactive points, README A6a.2)" begin
    # 4-point chain, point 4 inactive: its bonds are broken (bond damage 0) and the mask
    # gives row and column 4 the identity
    function chain(n, broken)
        coordinates = reshape(collect(0.0:(n - 1)), 1, n)
        nlist = neighbor_list(coordinates, 1.0)
        data = (coordinates = coordinates, volume = ones(n),
                nlist_ptr = collect(Int64, internal_element_ptr(nlist)),
                nlist = collect(flatview(nlist)),
                bond_damage = [(i in broken || j in broken) ? 0.0 : 1.0
                               for i in 1:n for j in nlist[i]])
        return data
    end
    el = TestSpring{1}()
    data = chain(4, (4,))
    pattern, K, ke = build_system(el, data, 4)
    compute_elements!(ke, pattern, el, data)
    assemble!(K, pattern, ke; active = [true, true, true, false])
    K4 = Matrix(to_sparse(K))
    data3 = chain(3, ())
    pattern3, K3, ke3 = build_system(el, data3, 3)
    compute_elements!(ke3, pattern3, el, data3)
    assemble!(K3, pattern3, ke3)
    @test K4[1:3, 1:3] == Matrix(to_sparse(K3))
    @test K4[4, :] == [0.0, 0.0, 0.0, 1.0]
    @test K4[:, 4] == [0.0, 0.0, 0.0, 1.0]
    # ghost rows (beyond n_own) of inactive points get no identity: the owner has it
    assemble!(K, pattern, ke; active = [true, true, true, false], n_own = 3)
    @test Matrix(to_sparse(K))[4, 4] == 0.0
    # pattern unchanged, active again: the plain assembly
    assemble!(K, pattern, ke; active = trues(4))
    @test Matrix(to_sparse(K))[4, 4] == 0.0 # bonds still broken
end
