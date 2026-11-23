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

# the former 8-loop rotation, kept as reference
function reference_fourth_order_rotation(R, C)
    n = size(C, 1)
    C_rot = zeros(n, n, n, n)
    for m in 1:n, nn in 1:n, o in 1:n, p in 1:n, i in 1:n, j in 1:n, k in 1:n, l in 1:n
        C_rot[m, nn, o, p] += C[i, j, k, l] * R[m, i] * R[nn, j] * R[o, k] * R[p, l]
    end
    return C_rot
end

function orthotropic_voigt(dof)
    dof == 2 && return [120.0 30.0 0.0; 30.0 40.0 0.0; 0.0 0.0 15.0]
    C = zeros(6, 6)
    C[1:3, 1:3] = [150.0 20.0 25.0; 20.0 60.0 30.0; 25.0 30.0 40.0]
    C[4, 4], C[5, 5], C[6, 6] = 10.0, 12.0, 14.0
    return C
end

@testset "fourth_order_rotation equals the 8-loop version" begin
    GZ = PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control
    for (dof, angles) in ((2, [33.0]), (3, [10.0, -20.0, 30.0]))
        C = copy(PeriLab.Helpers.get_fourth_order(orthotropic_voigt(dof), dof))
        R = PeriLab.Geometry.rotation_tensor(angles, dof)[1:dof, 1:dof]
        @test GZ.fourth_order_rotation(R, C) ≈ reference_fourth_order_rotation(R, C)
        @test GZ.fourth_order_rotation(R', C) ≈ reference_fourth_order_rotation(R', C)
    end
end

@testset "rotated zero energy stiffness equals the per-step rotation" begin
    GZ = PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control
    for dof in (2, 3)
        angle_sets = dof == 2 ? [[0.0], [30.0], [90.0], [-47.0]] :
                     [[0.0, 0.0, 0.0], [0.0, 0.0, 30.0], [90.0, 0.0, 0.0],
            [15.0, -40.0, 75.0]]
        nnodes = length(angle_sets)
        nodes = collect(1:nnodes)
        nv = size(orthotropic_voigt(dof), 1)
        CVoigt = zeros(nnodes, nv, nv)
        angles = zeros(nnodes, length(angle_sets[1]))
        Kinv = zeros(nnodes, dof, dof)
        for iID in nodes
            CVoigt[iID, :, :] = orthotropic_voigt(dof) .* iID
            angles[iID, :] = angle_sets[iID]
            K = rand(dof, dof)
            Kinv[iID, :, :] = K + K' + dof * [i == j for i in 1:dof, j in 1:dof]
        end
        z_ref = zeros(nnodes, dof, dof)
        for iID in nodes
            C = PeriLab.Helpers.get_fourth_order(CVoigt[iID, :, :], dof)
            R = PeriLab.Geometry.rotation_tensor(angles[iID, :], dof)[1:dof, 1:dof]
            GZ.global_zero_energy_mode_stiffness(iID, reference_fourth_order_rotation(R, C),
                                                 Kinv, z_ref)
        end
        C_rot = zeros(nnodes, nv, nv)
        GZ.rotate_voigt!(nodes, dof, CVoigt, angles, C_rot)
        z = zeros(nnodes, dof, dof)
        GZ.zero_energy_mode_stiffness!(nodes, C_rot, Kinv, z, Val(dof))
        @test isapprox(z, z_ref; rtol = 1e-12)
        @test z[1, :, :] ≈ z[1, :, :]'  # 0 deg: unrotated, symmetric Kinv and C
        for iID in nodes
            @test C_rot[iID, :, :] ≈ C_rot[iID, :, :]'
        end
    end
end
