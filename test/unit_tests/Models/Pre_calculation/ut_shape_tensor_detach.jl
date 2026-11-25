# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "Detach Damage" begin
    DM = PeriLab.Data_Manager
    Shape_Tensor = PeriLab.Solver_Manager.Model_Factory.Pre_Calculation.Shape_Tensor
    material(extra) = typed_block_material(merge(Dict("Material Model" => "Correspondence Elastic",
                                                      "Symmetry" => "isotropic",
                                                      "Bulk Modulus" => 1.0,
                                                      "Shear Modulus" => 1.0), extra))
    DM.initialize_data()
    DM.set_num_controller(3)
    DM.set_dof(2)
    damage = DM.create_node_scalar_field("Damage", Float64)[2]
    damage .= [0.2, 0.5, 0.9]
    Kinv = ones(3, 2, 2)

    # default 1: off, nothing changes
    m = material(Dict())
    @test m.base.detach_damage == 1.0
    DM.set_block_models(1, DM.BlockModels(material = m))
    Shape_Tensor.detach_nodes!(Kinv, [1, 2, 3], 1)
    @test Kinv == ones(3, 2, 2)

    # damage index >= value: K⁻¹ = 0; only the given nodes
    DM.set_block_models(1,
                        DM.BlockModels(material = material(Dict("Detach Damage" => 0.5))))
    Shape_Tensor.detach_nodes!(Kinv, [1, 2], 1)
    @test all(iszero, Kinv[2, :, :])
    @test Kinv[1, :, :] == ones(2, 2) && Kinv[3, :, :] == ones(2, 2)

    @test_throws PeriLab.PeriLabError material(Dict("Detach Damage" => 0.0))
    @test_throws PeriLab.PeriLabError material(Dict("Detach Damage" => 1.5))
end
