# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_Factory.jl")

using Test
include("../../../src/Support/data_manager.jl")
@testset "ut_valid_models" begin
    @test isnothing(FEM.valid_models(Dict()))
    @test isnothing(FEM.valid_models(Dict("Additive Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Damage Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Damage Model" => "a", "Additive Model" => "a", "Material Model" => "a Correspondence")))
    @test isnothing(FEM.valid_models(Dict("Thermal Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Material Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Material Model" => "a Correspondence")))
end

@testset "ut_init_FEM" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_dof(2)
    test_Data_manager.set_num_elements(3)
    test = FEM.init_FEM(test_Data_manager, Dict())
    @test isnothing(test)
    test_Data_manager.set_dof(1)
    test = FEM.init_FEM(test_Data_manager, Dict("Degree" => 1))
    @test isnothing(test)
    test = FEM.init_FEM(test_Data_manager, Dict("Degree" => 4))
    @test isnothing(test)
    test_Data_manager.set_dof(2)
    test_Data_manager = FEM.init_FEM(test_Data_manager, Dict("Degree" => 1, "Element Type" => "Lagrange"))
    @test "N Matrix" in test_Data_manager.get_all_field_keys()
    @test "B Matrix" in test_Data_manager.get_all_field_keys()
    N = test_Data_manager.get_field("N Matrix")
    @test size(N) == (4, 8, 2)
    @test N[1, 1:4, :] == [0.6220084679281462 0.0; 0.0 0.6220084679281462; 0.16666666666666663 0.0; 0.0 0.16666666666666663]
    B = test_Data_manager.get_field("B Matrix")
    @test size(B) == (4, 8, 3)
    @test B[3, 1:6, 3] == [-0.39433756729740643, -0.10566243270259354, -0.10566243270259354, 0.10566243270259354, 0.39433756729740643, -0.39433756729740643]
end
@testset "ut_get_polynomial_degree" begin

    @test isnothing(FEM.get_polynomial_degree(Dict(), 1))
    @test isnothing(FEM.get_polynomial_degree(Dict(), 2))
    @test isnothing(FEM.get_polynomial_degree(Dict(), 3))

    params = Dict("Degree" => 1)

    @test FEM.get_polynomial_degree(params, 2) == [1, 1]
    @test FEM.get_polynomial_degree(params, 3) == [1, 1, 1]

    params = Dict("Degree" => 2)
    @test FEM.get_polynomial_degree(params, 2) == [2, 2]
    @test FEM.get_polynomial_degree(params, 3) == [2, 2, 2]

    params = Dict("Degree" => 2.1)
    println()
    @test FEM.get_polynomial_degree(params, 2) == [2, 2]
    @test FEM.get_polynomial_degree(params, 3) == [2, 2, 2]

    params = Dict("Degree" => [2 3 1])
    @test isnothing(FEM.get_polynomial_degree(params, 2))
    @test FEM.get_polynomial_degree(params, 3) == [2, 3, 1]

    params = Dict("Degree" => [2.1 2])
    @test FEM.get_polynomial_degree(params, 2) == [2, 2]
    @test isnothing(FEM.get_polynomial_degree(params, 3))

end
