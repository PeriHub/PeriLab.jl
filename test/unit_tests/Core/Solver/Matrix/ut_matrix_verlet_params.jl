# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "ut_matrix_verlet_params" begin
    P = PeriLab.Solver_Manager.Matrix_Verlet.VerletMatrixParams
    p = typed_section(P, Dict{String,Any}())
    @test p.matrix_backend == "SparseArrays"
    @test p.gpu_weight == 1.0
    p = typed_section(P, Dict{String,Any}("Matrix Backend" => "Mixed", "GPU Weight" => 16))
    @test p.matrix_backend == "Mixed"
    @test p.gpu_weight == 16.0
    @test_throws PeriLab.PeriLabExceptions.PeriLabError typed_section(P,
                                                                      Dict{String,Any}("Matrix Backend" => "cpu"))
    reduction = Dict{String,Any}("Type" => "Static Condensation",
                                 "Reduction Blocks" => "2 3")
    @test typed_section(P, Dict{String,Any}("Model Reduction" => reduction)).model_reduction !==
          nothing
    @test_throws "Model Reduction requires Matrix Backend: SparseArrays" typed_section(P,
                                                                                       Dict{String,
                                                                                            Any}("Model Reduction" =>
                                                                                                     reduction,
                                                                                                 "Matrix Backend" => "CPU"))
end
