# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# One row of nodes, diffusivity calibrated to 1, foil |x| < 0.5 with initial concentration 1.
function read_stvd_results(filename)
    exo = ExodusDatabase(filename * ".e", "r")
    x = read_coordinates(exo)[1, :]
    n_steps = read_number_of_time_steps(exo)
    times = read_times(exo)
    # the first output is written before the initial condition is applied
    u = [read_values(exo, NodalVariable, step, 1, "Temperature") for step in 2:n_steps]
    close(exo)
    return x, times[2:end], u
end

second_moment(x, u) = sum(x .^ 2 .* u) / sum(u)

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    @testset "pure_diffusion" begin
        PeriLab.run("stvd_pure_diffusion.yaml"; silent = true)
        x, times, u = read_stvd_results("stvd_pure_diffusion")
        # the variance of the discrete solution grows exactly with 2 D t
        for i in eachindex(u)
            @test second_moment(x, u[i]) - second_moment(x, u[1])≈2 * (times[i] - times[1]) rtol=1e-4
            @test sum(u[i])≈sum(u[1]) rtol=1e-12
        end
        run_perilab("stvd_pure_diffusion", 1, true, folder_name)
    end
    @testset "linear" begin
        PeriLab.run("stvd_linear.yaml"; silent = true)
        x, times, u = read_stvd_results("stvd_linear")
        # front depth 0.5 * 0.5 = 0.25 -> core |x| < 0.25 is not reached yet
        core = abs.(x) .< 0.25
        @test all(u[end][core] .== 1.0)
        @test all(u[end][.!core] .< 1.0)
        @test sum(u[end])≈sum(u[1]) rtol=1e-12
        run_perilab("stvd_linear", 1, true, folder_name)
    end
    @testset "exponential" begin
        PeriLab.run("stvd_exponential.yaml"; silent = true)
        x, times, u = read_stvd_results("stvd_exponential")
        # front depth 0.5 * (1 - exp(-0.5)) = 0.197 -> core |x| < 0.303 is not reached yet
        core = abs.(x) .< 0.5 - 0.5 * (1 - exp(-0.5))
        @test all(u[end][core] .== 1.0)
        @test all(u[end][.!core] .< 1.0)
        @test sum(u[end])≈sum(u[1]) rtol=1e-12
        run_perilab("stvd_exponential", 1, true, folder_name)
    end
end
