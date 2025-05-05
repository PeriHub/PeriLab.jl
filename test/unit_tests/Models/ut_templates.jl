
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../../../src/Core/Data_manager.jl")
include("../../../src/Models/Degradation/Degradation_template/degradation_template.jl")
include("../../../src/Models/Additive/Additive_template/additive_template.jl")
include("../../../src/Models/Contact/Contact_template/contact_template.jl")
include("../../../src/Models/Damage/Damage_template/damage_template.jl")
include("../../../src/FEM/FEM_template/FEM_template.jl")
include("../../../src/Models/Material/Material_template/material_template.jl")
include("../../../src/Models/Material/Material_template/correspondence_template.jl")
include("../../../src/Models/Thermal/Thermal_template/thermal_template.jl")
include("../../../src/Models/Pre_calculation/Pre_calculation_template/pre_calculation_template.jl")

using Test
using TimerOutputs
#include("../../../src/PeriLab.jl")
#using .PeriLab
const to = TimerOutput()

test_data_manager = PeriLab.Data_manager
test_data_manager.initialize_data()
test_data_manager.set_num_controller(3)

@testset "ut_additive_template" begin
    @test Additive_template.additive_name() == "Additive Template"
    @test Additive_template.compute_model(test_data_manager,
                                          Vector{Int64}(1:3),
                                          Dict(),
                                          1,
                                          0.0,
                                          0.0) == test_data_manager
    @test Additive_template.init_model(test_data_manager, Vector{Int64}(1:3), Dict(), 1) ==
          test_data_manager
    @test Additive_template.fields_for_local_synchronization(test_data_manager, "") ==
          test_data_manager
end

@testset "ut_contact_template" begin
    @test Contact_template.contact_model_name() == "Contact Template"
    @test Contact_template.compute_model(test_data_manager,
                                         Vector{Int64}(1:3),
                                         Dict(),
                                         1,
                                         0.0,
                                         0.0) == test_data_manager
    @test Contact_template.init_contact_model(test_data_manager, Vector{Int64}(1:3), Dict(),
                                              1) ==
          test_data_manager
end

@testset "ut_degradation_template" begin
    @test Degradation_template.degradation_name() == "Degradation Template"
    @test Degradation_template.compute_model(test_data_manager,
                                             Vector{Int64}(1:3),
                                             Dict(),
                                             1,
                                             0.0,
                                             0.0) == test_data_manager
    @test Degradation_template.init_model(test_data_manager, Vector{Int64}(1:3), Dict(),
                                          1) ==
          test_data_manager
    @test Degradation_template.fields_for_local_synchronization(test_data_manager, "") ==
          test_data_manager
end

@testset "ut_damage_template" begin
    @test Damage_template.damage_name() == "Damage Template"
    @test Damage_template.compute_model(test_data_manager,
                                        Vector{Int64}(1:3),
                                        Dict(),
                                        1,
                                        0.0,
                                        0.0,
                                        to) == test_data_manager

    @test Damage_template.init_model(test_data_manager, Vector{Int64}(1:3), Dict(), 1) ==
          test_data_manager
    @test Damage_template.fields_for_local_synchronization(test_data_manager, "") ==
          test_data_manager
end

@testset "ut_FEM_template" begin
    @test FEM_template.element_name() == "element Template"
    @test FEM_template.compute_element(test_data_manager,
                                       Vector{Int64}(1:3),
                                       Dict(),
                                       0.0,
                                       0.0) == test_data_manager
    @test FEM_template.init_element(test_data_manager, Vector{Int64}(1:3), Dict(), [1]) ==
          test_data_manager
end

@testset "ut_material_template" begin
    test_data_manager = PeriLab.Data_manager
    @test !(Material_template.fe_support())
    @test Material_template.material_name() == "Material Template"
    @test Material_template.init_model(test_data_manager, Vector{Int64}(1:3), Dict()) ==
          test_data_manager
    @test Material_template.compute_model(test_data_manager,
                                          Vector{Int64}(1:3),
                                          Dict(),
                                          1,
                                          0.0,
                                          0.0,
                                          to) == test_data_manager
    @test Material_template.fields_for_local_synchronization(test_data_manager, "") ==
          test_data_manager
end

@testset "ut_thermal_template" begin
    @test Thermal_template.thermal_model_name() == "Thermal Template"
    @test Thermal_template.compute_model(test_data_manager,
                                         Vector{Int64}(1:3),
                                         Dict(),
                                         1,
                                         0.0,
                                         0.0) == test_data_manager
    @test Thermal_template.init_model(test_data_manager, Vector{Int64}(1:3), Dict()) ==
          test_data_manager
end

@testset "ut_correspondence_template" begin
    test_data_manager = PeriLab.Data_manager
    @test !(Correspondence_template.fe_support())
    @test Correspondence_template.correspondence_name() == "Correspondence Template"

    @test Correspondence_template.init_model(test_data_manager,
                                             Vector{Int64}(1:3),
                                             Dict(),
                                             1) == test_data_manager

    dat,
    vec = Correspondence_template.compute_stresses(test_data_manager,
                                                   1,
                                                   2,
                                                   Dict(),
                                                   0.0,
                                                   0.0,
                                                   view([1, 2], :, :, :),
                                                   view([1, 0], :, :, :),
                                                   view([-1, 2.2], :, :, :))
    @test dat == test_data_manager
    @test vec[1] == -1
    @test vec[2] == 2.2

    dat,
    vec = Correspondence_template.compute_stresses(test_data_manager,
                                                   1,
                                                   2,
                                                   Dict(),
                                                   0.0,
                                                   0.0,
                                                   view([1, 2], :, :, :),
                                                   view([1, 0], :, :, :),
                                                   view([-1, 2.2], :, :, :),
                                                   (1, 1))
    @test dat == test_data_manager
    @test vec[1] == -1
    @test vec[2] == 2.2

    vec,
    dat = Correspondence_template.compute_stresses(test_data_manager,
                                                   2,
                                                   Dict(),
                                                   0.0,
                                                   0.0,
                                                   [1.0, 2.0],
                                                   [1.0, 0.0],
                                                   [-1.0, 2.2])
    @test dat == test_data_manager
    @test vec[1] == -1
    @test vec[2] == 2.2
    @test Correspondence_template.fields_for_local_synchronization(test_data_manager, "") ==
          test_data_manager
end

@testset "ut_pre_calculation_template" begin
    @test Pre_calculation_template.pre_calculation_name() == "pre_calculation Template"
    @test Pre_calculation_template.compute_model(test_data_manager,
                                                 Vector{Int64}(1:3),
                                                 Dict(),
                                                 1) == test_data_manager
    @test Pre_calculation_template.init_model(test_data_manager,
                                              Vector{Int64}(1:3),
                                              Dict(),
                                              1) == test_data_manager
    @test Pre_calculation_template.fields_for_local_synchronization(test_data_manager,
                                                                    "") == test_data_manager
end
