
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../../../src/Support/data_manager.jl")
include("../../../src/Physics/Corrosion/Corrosion_template/corrosion_template.jl")
include("../../../src/Physics/Additive/Additive_template/additive_template.jl")
include("../../../src/Physics/Damage/Damage_template/damage_template.jl")
include("../../../src/Physics/Material/Material_template/material_template.jl")
include("../../../src/Physics/Material/Material_template/correspondence_template.jl")
include("../../../src/Physics/Thermal/Thermal_template/thermal_template.jl")
include("../../../src/Physics/Pre_calculation/Pre_calculation_template/pre_calculation_template.jl")

using Test
using TimerOutputs
#include("../../../src/PeriLab.jl")
#using .PeriLab
const to = TimerOutput()

test_Data_manager = PeriLab.Data_manager
test_Data_manager.set_num_controller(3)

@testset "ut_additive_template" begin
    @test Additive_template.additive_name() == "Additive Template"
    @test Additive_template.compute_additive_model(test_Data_manager, Vector{Int64}(1:3), Dict(), 0.0, 0.0) == test_Data_manager
    @test Additive_template.init_additive_model(test_Data_manager, Vector{Int64}(1:3), Dict(), 1) == test_Data_manager
end

@testset "ut_additive_template" begin
    @test Corrosion_template.corrosion_name() == "Corrosion Template"
    @test Corrosion_template.compute_corrosion_model(test_Data_manager, Vector{Int64}(1:3), Dict(), 0.0, 0.0) == test_Data_manager
    @test Corrosion_template.init_corrosion_model(test_Data_manager, Vector{Int64}(1:3), Dict(), 1) == test_Data_manager
end

@testset "ut_damage_template" begin
    @test Damage_template.damage_name() == "Damage Template"
    @test Damage_template.compute_damage(test_Data_manager, Vector{Int64}(1:3), Dict(), 1, 0.0, 0.0) == test_Data_manager
    @test Damage_template.compute_damage_pre_calculation(test_Data_manager, Vector{Int64}(1:3), 1, "dummy for function", 0.0, 0.0) == test_Data_manager
    # @test Damage_template.init_damage_model(test_Data_manager, Vector{Int64}(1:3), Dict(), 1) == test_Data_manager
end

@testset "ut_material_template" begin
    test_Data_manager = PeriLab.Data_manager
    @test !(Material_template.fe_support())
    @test Material_template.material_name() == "Material Template"
    @test Material_template.init_material_model(test_Data_manager, Vector{Int64}(1:3), Dict()) == test_Data_manager
    @test Material_template.compute_forces(test_Data_manager, Vector{Int64}(1:3), Dict(), 0.0, 0.0, to) == test_Data_manager
end

@testset "ut_correspondence_template" begin
    test_Data_manager = PeriLab.Data_manager
    @test !(Correspondence_template.fe_support())
    @test Correspondence_template.correspondence_name() == "Correspondence Template"
    @test Correspondence_template.init_material_model(test_Data_manager, Vector{Int64}(1:3), Dict()) == test_Data_manager

    dat, vec = Correspondence_template.compute_stresses(test_Data_manager, 1, 2, Dict(), 0.0, 0.0, view([1, 2], :, :, :), view([1, 0], :, :, :), view([-1, 2.2], :, :, :))
    @test dat == test_Data_manager
    @test vec[1] == -1
    @test vec[2] == 2.2
end

@testset "ut_thermal_template" begin
    @test Thermal_template.thermal_model_name() == "Thermal Template"
    @test Thermal_template.compute_thermal_model(test_Data_manager, Vector{Int64}(1:3), Dict(), 0.0, 0.0) == test_Data_manager
    @test Thermal_template.init_thermal_model(test_Data_manager, Vector{Int64}(1:3), Dict()) == test_Data_manager
end

@testset "ut_pre_calculation_template" begin
    @test Pre_calculation_template.pre_calculation_name() == "pre_calculation Template"
    @test Pre_calculation_template.pre_calculation(test_Data_manager, Vector{Int64}(1:3), Dict(), 0.0, 0.0) == test_Data_manager
end
