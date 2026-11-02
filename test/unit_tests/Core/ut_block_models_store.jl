# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "empty block models" begin
    PeriLab.Data_Manager.initialize_data()
    m = PeriLab.Data_Manager.get_block_models(1)
    @test m.material === nothing && m.damage === nothing && m.thermal === nothing
    @test m.additive === nothing && m.degradation === nothing
    @test isempty(m.pre_calculation)
    MF = PeriLab.Solver_Manager.Model_Factory
    for name in ("Material Model", "Damage Model", "Thermal Model", "Additive Model",
                 "Degradation Model", "Pre Calculation Model")
        @test !MF.has_block_model(1, name)
    end
    @test MF.local_damping_symmetry(1) == "3D"
    @test MF.block_local_damping(1) === nothing
    @test MF.block_thermal_conductivity(1) === nothing
end

@testset "check_dependencies updates block models" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(3)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                                  "Shear Modulus" => 1.0))
    th = typed_model(:thermal, Dict("Thermal Model" => "Heat Transfer",
                                    "Heat Transfer Coefficient" => 1.0,
                                    "Environmental Temperature" => 30);
                     name_key = "Thermal Model")
    PeriLab.Data_Manager.set_block_models(1,
                                          PeriLab.Data_Manager.BlockModels(material = m,
                                                                           thermal = th,
                                                                           pre_calculation = ["Shape Tensor"]))
    PeriLab.Solver_Manager.Model_Factory.Pre_Calculation.check_dependencies(Dict(1 => [1, 2]))
    b = PeriLab.Data_Manager.get_block_models(1)
    @test b.pre_calculation == ["Deformed Bond Geometry", "Shape Tensor", "Deformation Gradient"]
    @test b.material === m && b.thermal === th             # other parts kept
end
