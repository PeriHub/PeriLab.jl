# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using PeriLab
#using Test

# from Peridigm

nnodes = 5
dof = 2

PeriLab.Data_Manager.initialize_data()
PeriLab.Data_Manager.set_num_controller(5)
PeriLab.Data_Manager.set_dof(2)
blocks = PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64)
horizon = PeriLab.Data_Manager.create_constant_node_scalar_field("Horizon", Float64)
coor = PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates", Float64, 2)
density = PeriLab.Data_Manager.create_constant_node_scalar_field("Density", Float64)
volume = PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64)
length_nlist = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                      Int64)
length_nlist .= 4

nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
undeformed_bond = PeriLab.Data_Manager.create_constant_bond_vector_state("Bond Geometry",
                                                                         Float64,
                                                                         dof)
undeformed_bond_length = PeriLab.Data_Manager.create_constant_bond_scalar_state("Bond Length",
                                                                                Float64)
heat_capacity = PeriLab.Data_Manager.create_constant_node_scalar_field("Specific Heat Capacity",
                                                                       Float64;
                                                                       default_value = 18000)
nlist[1] = [2, 3, 4, 5]
nlist[2] = [1, 3, 4, 5]
nlist[3] = [1, 2, 4, 5]
nlist[4] = [1, 2, 3, 5]
nlist[5] = [1, 2, 3, 4]

coor[1, 1] = 0;
coor[1, 2] = 0;
coor[2, 1] = 0.5;
coor[2, 2] = 0.5;
coor[3, 1] = 1;
coor[3, 2] = 0;
coor[4, 1] = 0;
coor[4, 2] = 1;
coor[5, 1] = 1;
coor[5, 2] = 1;

volume = [0.5, 0.5, 0.5, 0.5, 0.5]
density = [1e-6, 1e-6, 3e-6, 3e-6, 1e-6]
horizon = [3.1, 3.1, 3.1, 3.1, 3.1]

PeriLab.Geometry.bond_geometry!(undeformed_bond,
                                undeformed_bond_length,
                                Vector(1:nnodes),
                                nlist,
                                coor)

blocks = ["1", "2"]
blocks = PeriLab.Data_Manager.set_block_name_list(blocks)
@testset "ut_mechanical_critical_time_step" begin
    t = PeriLab.Solver_Manager.Model_Factory.compute_mechanical_critical_time_step(Vector{Int64}(1:nnodes),
                                                                                   Float64(140.0))
    @test t == 1.4142135623730952e25 # not sure if this is right :D
end
# from Peridigm
@testset "ut_thermodynamic_crititical_time_step" begin
    t = PeriLab.Solver_Manager.Model_Factory.compute_thermodynamic_critical_time_step(Vector{Int64}(1:nnodes),
                                                                                      Float64(0.12))
    @test t == 1e25
end

@testset "ut_read_properties" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_block_name_list(["block_1"])
    PeriLab.Data_Manager.set_block_id_list([1])
    input = typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                "Density" => 1.0,
                                                                "Horizon" => 1.0,
                                                                "Thermal Model" => "therm")),
                             "Models" => Dict("Thermal Models" => Dict("therm" => Dict("Thermal Model" => "Heat Transfer",
                                                                                       "Heat Transfer Coefficient" => 1.0,
                                                                                       "Environmental Temperature" => 30)))))
    PeriLab.Solver_Manager.Model_Factory.read_properties(input, false)
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 1) === input.thermals["therm"]
end

@testset "ut_add_model" begin
    PeriLab.Data_Manager.initialize_data()
    @test_logs (:error,
                "Model Test is not specified and cannot be included.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Model_Factory.add_model("Test")
    end
end

@testset "ut_test_timestep" begin
    @test PeriLab.Solver_Manager.Model_Factory.test_timestep(1.0, 2.0) == 1
    @test PeriLab.Solver_Manager.Model_Factory.test_timestep(2.0, 1.1) == 1.1
    @test PeriLab.Solver_Manager.Model_Factory.test_timestep(2.0, 2.0) == 2
end

@testset "ut_get_cs_denominator" begin
    volume = Float64[1, 2, 3]
    undeformed_bond = [1.0, 2, 3]
    @test PeriLab.Solver_Manager.Model_Factory.get_cs_denominator(volume,
                                                                  undeformed_bond) == 3
    undeformed_bond = [2.0, 4, 6]
    @test PeriLab.Solver_Manager.Model_Factory.get_cs_denominator(volume,
                                                                  undeformed_bond) == 1.5
    undeformed_bond = [1.0, 0.5, 2]
    @test PeriLab.Solver_Manager.Model_Factory.get_cs_denominator(volume,
                                                                  undeformed_bond) == 6.5
end

@testset "read_properties builds typed block materials" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    input = typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                 "Density" => 1.0,
                                                                 "Horizon" => 1.0,
                                                                 "Material Model" => "Mat")),
                             "Models" => Dict("Material Models" => Dict("Mat" => Dict("Material Model" => "PD Solid Elastic",
                                                                                       "Symmetry" => "isotropic plane strain",
                                                                                       "Bulk Modulus" => 10.0,
                                                                                       "Shear Modulus" => 10.0)))))
    PeriLab.Data_Manager.set_block_name_list(["block_1"])
    PeriLab.Data_Manager.set_block_id_list([1])
    PeriLab.Solver_Manager.Model_Factory.read_properties(input, true)
    m = PeriLab.Data_Manager.get_block_material(1)
    @test m isa PeriLab.Solver_Manager.Model_Factory.Material.BlockMaterial
    @test m.symmetry == "plane strain"
    @test m.moduli.youngs_modulus == 22.5

end


@testset "read_properties aborts on an undefined material name" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_block_name_list(["block_1"])
    PeriLab.Data_Manager.set_block_id_list([1])
    models = Dict("Material Models" => Dict("Mat" => Dict("Material Model" => "PD Solid Elastic",
                                                          "Symmetry" => "isotropic plane strain",
                                                          "Bulk Modulus" => 10.0,
                                                          "Shear Modulus" => 10.0)))
    input = typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                 "Density" => 1.0,
                                                                 "Horizon" => 1.0,
                                                                 "Material Model" => "Steeel")),
                             "Models" => models))
    @test_logs (:error,
                "Material Model model with name Steeel is defined in blocks, but missing in the Material Models definition.") @test_throws PeriLab.PeriLabError PeriLab.Solver_Manager.Model_Factory.read_properties(input,
                                                                                                                                                                                                                     true)
    # without a material model the name is not checked
    PeriLab.Solver_Manager.Model_Factory.read_properties(input, false)
    @test PeriLab.Data_Manager.get_block_material(1) === nothing
end

@testset "local damping symmetry of a block without material" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(2)
    @test PeriLab.Solver_Manager.Model_Factory.local_damping_symmetry(1) == "3D"
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Symmetry" => "isotropic plane strain",
                                  "Bulk Modulus" => 10.0, "Shear Modulus" => 10.0))
    PeriLab.Data_Manager.set_block_material(1, m)
    @test PeriLab.Solver_Manager.Model_Factory.local_damping_symmetry(1) == "plane strain"
end

@testset "read_properties builds block damages" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_block_name_list(["block_1", "block_2"])
    PeriLab.Data_Manager.set_block_id_list([1, 2])
    blocks = Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0, "Horizon" => 1.0,
                                    "Damage Model" => "Dam"),
                  "block_2" => Dict("Block ID" => 2, "Density" => 1.0, "Horizon" => 1.0))
    models = Dict("Damage Models" => Dict("Dam" => Dict("Damage Model" => "Critical Stretch",
                                                        "Critical Value" => 0.1)))
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    PeriLab.Solver_Manager.Model_Factory.read_properties(input, false)
    d = PeriLab.Data_Manager.get_block_damage(1)
    @test d isa PeriLab.Solver_Manager.Model_Factory.Damage.BlockDamage
    @test PeriLab.Data_Manager.get_block_damage(2) === nothing

    blocks["block_1"]["Damage Model"] = "Dmg"
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    @test_logs (:error,
                "Damage Model model with name Dmg is defined in blocks, but missing in the Damage Models definition.") match_mode=:any @test_throws PeriLab.PeriLabError PeriLab.Solver_Manager.Model_Factory.read_properties(input,
                                                                                                                                                                                                                                false)
    input = typed_input(Dict("Blocks" => blocks, "Models" => Dict()))
    @test_logs (:error,
                "Damage Model is defined in blocks, but no Damage Models definition block exists") match_mode=:any @test_throws PeriLab.PeriLabError PeriLab.Solver_Manager.Model_Factory.read_properties(input,
                                                                                                                                                                                                                 false)
end

@testset "local damping reads the block damage" begin
    PeriLab.Data_Manager.initialize_data()
    MF = PeriLab.Solver_Manager.Model_Factory
    @test !MF.has_block_model(1, "Damage Model")
    damping = Dict("Representative Young's modulus" => 1.0, "Damping coefficient" => 0.5)
    damages = Dict("Dam" => Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                                 "Local Damping" => damping))
    blocks = Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0, "Horizon" => 1.0,
                                    "Damage Model" => "Dam"))
    input = typed_input(Dict("Blocks" => blocks, "Models" => Dict("Damage Models" => damages)))
    d = MF.Damage.block_damage(input.damages["Dam"])
    PeriLab.Data_Manager.set_block_damage(1, d)
    @test MF.has_block_model(1, "Damage Model")
    @test MF.block_model_parameters(1, "Damage Model") === d
    @test MF.block_local_damping(1).damping_coefficient == 0.5
    @test MF.block_local_damping(2) === nothing
end
