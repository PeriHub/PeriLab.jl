# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

@testset "control_name" begin
    @test PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control.control_name() ==
          "Global"
end

@testset "rotate_fourth_order_tensor_interface test" begin
    testval = zeros(3, 3, 3, 3)
    @test PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control.rotate_fourth_order_tensor(zeros(3),
                                                                                                           zeros(3,
                                                                                                                 3,
                                                                                                                 3,
                                                                                                                 3),
                                                                                                           3,
                                                                                                           true) ==
          testval
    @test PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control.rotate_fourth_order_tensor(zeros(3),
                                                                                                           zeros(3,
                                                                                                                 3,
                                                                                                                 3,
                                                                                                                 3),
                                                                                                           3,
                                                                                                           false) ==
          testval
    testval = zeros(2, 2, 2, 2)
    @test PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control.rotate_fourth_order_tensor(zeros(1),
                                                                                                           zeros(2,
                                                                                                                 2,
                                                                                                                 2,
                                                                                                                 2),
                                                                                                           2,
                                                                                                           true) ==
          testval
    @test PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control.rotate_fourth_order_tensor(zeros(1),
                                                                                                           zeros(2,
                                                                                                                 2,
                                                                                                                 2,
                                                                                                                 2),
                                                                                                           2,
                                                                                                           false) ==
          testval
    @test PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control.rotate_fourth_order_tensor(zeros(3),
                                                                                                           zeros(2,
                                                                                                                 2,
                                                                                                                 2,
                                                                                                                 2),
                                                                                                           3,
                                                                                                           true) ==
          testval
    @test PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control.rotate_fourth_order_tensor(zeros(3),
                                                                                                           zeros(2,
                                                                                                                 2,
                                                                                                                 2,
                                                                                                                 2),
                                                                                                           3,
                                                                                                           false) ==
          testval
end

@testset "rotate_fourth_order_tensor leaves the input unchanged" begin
    GZ = PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control
    C = zeros(3, 3, 3, 3)
    C[1, 1, 1, 1] = 3.0
    C[2, 2, 2, 2] = 1.0
    C[1, 2, 1, 2] = C[1, 2, 2, 1] = C[2, 1, 1, 2] = C[2, 1, 2, 1] = 0.5
    C_in = copy(C)
    C_rot = GZ.rotate_fourth_order_tensor([0.0, 0.0, 30.0], C, 3, false)
    @test C == C_in
    @test GZ.rotate_fourth_order_tensor([0.0, 0.0, 30.0], C_rot, 3, true) ≈ C
    @test GZ.rotate_fourth_order_tensor([0.0, 0.0, 90.0], C, 3, false)[2, 2, 2, 2] ≈ 3.0
end
