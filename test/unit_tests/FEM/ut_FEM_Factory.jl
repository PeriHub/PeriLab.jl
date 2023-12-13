# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_Factory.jl")

using Test
include("../../../src/Support/data_manager.jl")


test_Data_manager = Data_manager
dof = 2
test_Data_manager.set_dof(dof)
test_Data_manager.set_num_elements(2)
test_Data_manager.set_num_controller(6)
test_Data_manager.create_node_field("Force Density", Float64, dof)
test_Data_manager.create_node_field("Displacements", Float64, dof)
test_Data_manager.create_constant_node_field("Coordinates", Float64, dof)
params = Dict("FEM" => Dict("FE_1" => Dict("Degree" => 1, "Element Type" => "Lagrange", "Material Model" => "Elastic Model")),
    "Material Models" => Dict("Elastic Model" => Dict("Material Model" => "Correspondence Elastic", "Symmetry" => "isotropic plane strain", "Young's Modulus" => 2.5e+3, "Poisson's Ratio" => 0.33, "Shear Modulus" => 2.0e3)))


topology = test_Data_manager.create_constant_free_size_field("FE Element Topology", Int64, (2, 4))
topology[1, 1] = 1
topology[1, 2] = 2
topology[1, 3] = 3
topology[1, 4] = 4
topology[2, 1] = 3
topology[2, 2] = 5
topology[2, 3] = 4
topology[2, 4] = 6
test_Data_manager = FEM.init_FEM(test_Data_manager, params["FEM"]["FE_1"])
elements = Vector{Int64}([1, 2])
#test_Data_manager = FEM.eval(test_Data_manager, elements, params, "FE_1", 0.0, 1.0e-6)
@testset "ut_valid_models" begin
    @test isnothing(FEM.valid_models(Dict()))
    @test isnothing(FEM.valid_models(Dict("Additive Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Damage Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Damage Model" => "a", "Additive Model" => "a", "Material Model" => "a Correspondence")))
    @test isnothing(FEM.valid_models(Dict("Thermal Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Material Model" => "a")))
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

    @test "Element StrainN" in test_Data_manager.get_all_field_keys()
    @test "Element StressN" in test_Data_manager.get_all_field_keys()
    @test "Element StrainNP1" in test_Data_manager.get_all_field_keys()
    @test "Element StressNP1" in test_Data_manager.get_all_field_keys()
    @test "Element Strain Increment" in test_Data_manager.get_all_field_keys()

end



