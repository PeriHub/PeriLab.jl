# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


include("../../../../../src/Physics/Material/Material_Models/Zero_Energy_Control/global_control.jl")
using Test
using .global_zero_energy_control


@testset "control_name" begin
    @test global_zero_energy_control.control_name() == "Global"
end

@testset "rotate_fourth_order_tensor_interface test" begin
    testval = zeros(3, 3, 3, 3)
    @test global_zero_energy_control.rotate_fourth_order_tensor(zeros(3), zeros(3, 3, 3, 3), 3, true) == testval
    @test global_zero_energy_control.rotate_fourth_order_tensor(zeros(3), zeros(3, 3, 3, 3), 3, false) == testval
    testval = zeros(2, 2, 2, 2)
    @test global_zero_energy_control.rotate_fourth_order_tensor(zeros(1), zeros(2, 2, 2, 2), 2, true) == testval
    @test global_zero_energy_control.rotate_fourth_order_tensor(zeros(1), zeros(2, 2, 2, 2), 2, false) == testval
    @test global_zero_energy_control.rotate_fourth_order_tensor(zeros(3), zeros(2, 2, 2, 2), 2, true) == testval
    @test global_zero_energy_control.rotate_fourth_order_tensor(zeros(3), zeros(2, 2, 2, 2), 2, false) == testval
end