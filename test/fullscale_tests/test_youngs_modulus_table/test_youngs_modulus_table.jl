# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Young's Modulus given as a table over the Crystallinity field (E_over_X.txt,
# E = 2100 + 5110 X MPa). Box 40 x 4 x 2 mm, E_eff = (F/A) / strain between x = 15 and
# x = 25 with F the reaction force at the clamp.

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    E(X) = 2100 + 5110 * X
    dx = 0.5
    x = [(i - 0.5) * dx for i in 1:80, _ in 1:8, _ in 1:4]
    y = [(j - 0.5) * dx for _ in 1:80, j in 1:8, _ in 1:4]
    z = [(k - 0.5) * dx for _ in 1:80, _ in 1:8, k in 1:4]
    function write_mesh(crystallinity, block)
        open("box_mesh.txt", "w") do io
            println(io, "header: x y z block_id volume Crystallinity")
            for n in eachindex(x)
                println(io, "$(x[n]) $(y[n]) $(z[n]) $(block(x[n])) $(dx^3) ",
                        crystallinity(x[n]))
            end
        end
    end
    function run_box(deck = "box.yaml")
        @test PeriLab.run(deck; silent = true) == 0
        exo = ExodusDatabase("box.e", "r")
        n = read_number_of_time_steps(exo)
        u = read_values(exo, NodalVariable, n, 1, "Displacementsx")
        F = sum(read_values(exo, NodalVariable, n, 1, "Forcesx")[vec(x .< 1.5)])
        # the first output is written before the material is initialized
        E_field = [read_values(exo, NodalVariable, step, 1, "Young's Modulus")
                   for step in 2:n]
        close(exo)
        rm("box.e")
        return u, F, E_field
    end
    # mean strain between the cross sections x0 and x1 (node centers)
    function strain(u, x0, x1)
        mean_u(xs) = sum(u[vec(abs.(x .- xs) .< 1e-6)]) / 32
        return (mean_u(x1) - mean_u(x0)) / (x1 - x0)
    end

    @testset "homogeneous box, X = 0.18" begin
        write_mesh(_ -> 0.18, _ -> 1)
        u, F, E_field = run_box()
        @test E_field[end] ≈ fill(E(0.18), length(x))
        @test F / (4 * 2) / strain(u, 15.25, 24.75)≈E(0.18) rtol=0.02
    end

    @testset "two blocks, X = 0 | X = 0.5" begin
        write_mesh(xi -> xi < 20 ? 0.0 : 0.5, xi -> xi < 20 ? 1 : 2)
        u, _, _ = run_box()
        # away from the clamp, the load and the interface
        @test strain(u, 8.25, 15.25) / strain(u, 24.75, 31.75)≈E(0.5) / E(0.0) rtol=0.02
    end

    @testset "no compounding over solver steps" begin
        write_mesh(xi -> xi < 20 ? 0.0 : 0.5, xi -> xi < 20 ? 1 : 2)
        _, _, E_field = run_box("box_two_steps.yaml")
        expected = [xi < 20 ? E(0.0) : E(0.5) for xi in vec(x)]
        @test length(E_field) >= 2
        @test all(E_step ≈ expected for E_step in E_field)
    end
    rm("box_mesh.txt")
end
