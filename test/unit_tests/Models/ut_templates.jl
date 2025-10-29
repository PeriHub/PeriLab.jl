
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../../../src/Core/Data_Manager.jl")
include("../../../src/Models/Degradation/Degradation_template/degradation_template.jl")
include("../../../src/Models/Additive/Additive_template/additive_template.jl")
include("../../../src/Models/Contact/Contact_template/contact_template.jl")
include("../../../src/Models/Damage/Damage_template/damage_template.jl")
include("../../../src/FEM/FEM_template/FEM_template.jl")
include("../../../src/Models/Material/Material_Models/Material_template/material_template.jl")
include("../../../src/Models/Material/Material_Models/Material_template/correspondence_template.jl")
include("../../../src/Models/Thermal/Thermal_template/thermal_template.jl")
include("../../../src/Models/Pre_calculation/Pre_calculation_template/pre_calculation_template.jl")

using Test
#using PeriLab

test_data_manager = PeriLab.Data_Manager
test_data_manager.initialize_data()
test_data_manager.set_num_controller(3)

@testset "ut_additive_template" begin
    @test Additive_template.additive_name() == "Additive Template"
    Additive_template.compute_model(Vector{Int64}(1:3),
                                    Dict(),
                                    1,
                                    0.0,
                                    0.0)
    @test Additive_template.init_model(Vector{Int64}(1:3), Dict(), 1) ==
          test_data_manager
    @test Additive_template.fields_for_local_synchronization("") ==
          test_data_manager
end

@testset "ut_contact_template" begin
    @test Contact_template.contact_model_name() == "Contact Template"
    Contact_template.compute_model(Vector{Int64}(1:3),
                                   Dict(),
                                   1,
                                   0.0,
                                   0.0)
    @test Contact_template.init_contact_model(Vector{Int64}(1:3), Dict(),
                                              1) ==
          test_data_manager
end

@testset "ut_degradation_template" begin
    @test Degradation_template.degradation_name() == "Degradation Template"
    Degradation_template.compute_model(Vector{Int64}(1:3),
                                       Dict(),
                                       1,
                                       0.0,
                                       0.0)
    @test Degradation_template.init_model(Vector{Int64}(1:3), Dict(),
                                          1) ==
          test_data_manager
    @test Degradation_template.fields_for_local_synchronization("") ==
          test_data_manager
end

@testset "ut_damage_template" begin
    @test Damage_template.damage_name() == "Damage Template"
    Damage_template.compute_model(Vector{Int64}(1:3),
                                  Dict(),
                                  1,
                                  0.0,
                                  0.0)

    @test Damage_template.init_model(Vector{Int64}(1:3), Dict(), 1) ==
          test_data_manager
    @test Damage_template.fields_for_local_synchronization("") ==
          test_data_manager
end

@testset "ut_FEM_template" begin
    @test FEM_template.element_name() == "element Template"
    FEM_template.compute_element(Vector{Int64}(1:3),
                                 Dict(),
                                 0.0,
                                 0.0)
    @test FEM_template.init_element(Vector{Int64}(1:3), Dict(), [1]) ==
          test_data_manager
end

@testset "ut_material_template" begin
    test_data_manager = PeriLab.Data_Manager
    @test !(Material_template.fe_support())
    @test Material_template.material_name() == "Material Template"
    @test Material_template.init_model(Vector{Int64}(1:3), Dict()) ==
          test_data_manager
    Material_template.compute_model(Vector{Int64}(1:3),
                                    Dict(),
                                    1,
                                    0.0,
                                    0.0)
    @test Material_template.fields_for_local_synchronization("") ==
          test_data_manager
end

@testset "ut_thermal_template" begin
    @test Thermal_template.thermal_model_name() == "Thermal Template"
    Thermal_template.compute_model(Vector{Int64}(1:3),
                                   Dict(),
                                   1,
                                   0.0,
                                   0.0)
    @test Thermal_template.init_model(Vector{Int64}(1:3), Dict()) ==
          test_data_manager
end

@testset "ut_correspondence_template" begin
    test_data_manager = PeriLab.Data_Manager
    @test !(Correspondence_template.fe_support())
    @test Correspondence_template.correspondence_name() == "Correspondence Template"

    Correspondence_template.init_model(Vector{Int64}(1:3),
                                       Dict(),
                                       1)

    vec = Correspondence_template.compute_stresses(1,
                                                   2,
                                                   Dict(),
                                                   0.0,
                                                   0.0,
                                                   view([1.0, 2.0], :, :, :),
                                                   view([1, 0.0], :, :, :),
                                                   view([-1, 2.2], :, :, :))
    @test vec[1] == -1
    @test vec[2] == 2.2

    vec = Correspondence_template.compute_stresses(1,
                                                   2,
                                                   Dict(),
                                                   0.0,
                                                   0.0,
                                                   view([1.0, 2.0], :, :, :),
                                                   view([1.0, 0.0], :, :, :),
                                                   view([-1, 2.2], :, :, :),
                                                   (1, 1))
    @test vec[1] == -1
    @test vec[2] == 2.2

    vec = Correspondence_template.compute_stresses(2,
                                                   Dict(),
                                                   0.0,
                                                   0.0,
                                                   [1.0, 2.0],
                                                   [1.0, 0.0],
                                                   [-1.0, 2.2])
    @test vec[1] == -1
    @test vec[2] == 2.2
    Correspondence_template.fields_for_local_synchronization("")
end

@testset "ut_pre_calculation_template" begin
    @test Pre_calculation_template.pre_calculation_name() == "pre_calculation Template"
    Pre_calculation_template.compute_model(Vector{Int64}(1:3),
                                           Dict(),
                                           1)
    Pre_calculation_template.init_model(Vector{Int64}(1:3),
                                        Dict(),
                                        1)
    Pre_calculation_template.fields_for_local_synchronization("")
end
