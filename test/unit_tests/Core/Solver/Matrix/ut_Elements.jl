# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "ut_Elements" begin
    @test element_type("Test Spring") === nothing
    register_element("Test Spring", TestSpring{1})
    @test element_type("Test Spring") === TestSpring{1}
    @test "Test Spring" in registered_elements()
    delete!(PeriLab.Solver_Manager.PD_Matrix.Elements.REGISTRY, "Test Spring")
end

@testset "ut_Backend" begin
    comm = MPI.COMM_WORLD
    @test select_backend("CPU", comm) == (Array, false, -1)
    @test_throws ArgumentError select_backend("cpu", comm)
    PeriLab.Solver_Manager.PD_Matrix.Backend.load_cuda!()
    if GPU_HOOK[] === nothing
        @test_throws "Matrix Backend GPU requires CUDA.jl (using CUDA) and a functional GPU" select_backend("GPU",
                                                                                                            comm)
        @test select_backend("Mixed", comm) == (Array, false, -1)
    end
end
