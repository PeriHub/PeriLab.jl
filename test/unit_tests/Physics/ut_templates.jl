
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Support/data_manager.jl")
include("../../../src/Physics/Additive/Additive_template/additive_template.jl")
include("../../../src/Physics/Damage/Damage_template/damage_template.jl")
include("../../../src/Physics/Material/Material_template/material_template.jl")
include("../../../src/Physics/Thermal/Thermal_template/thermal_template.jl")

using Test

testDatamanager = Data_manager
testDatamanager.set_nmasters(3)

@testset "ut_additive_template" begin
    @test Additive_template.additive_name() == "Additive Template"
    @test Additive_template.compute_additive(testDatamanager, Vector{Int64}(1:3), Dict(), 0.0, 0.0) == testDatamanager

end

@testset "ut_damage_template" begin
    @test Damage_template.damage_name() == "Damage Template"
    @test Damage_template.compute_damage(testDatamanager, Vector{Int64}(1:3), Dict(), 1, 0.0, 0.0) == testDatamanager
    @test Damage_template.compute_damage_pre_calculation(testDatamanager, Vector{Int64}(1:3), 1, "dummy for function", 0.0, 0.0) == testDatamanager
end

@testset "ut_material_template" begin
    @test Material_template.material_name() == "Material Template"
    @test Material_template.compute_force(testDatamanager, Vector{Int64}(1:3), Dict(), 0.0, 0.0) == testDatamanager
end

@testset "ut_thermal_template" begin
    @test Thermal_template.thermal_model_name() == "Thermal Template"
    @test Thermal_template.compute_thermal_model(testDatamanager, Vector{Int64}(1:3), Dict(), 0.0, 0.0) == testDatamanager
end
