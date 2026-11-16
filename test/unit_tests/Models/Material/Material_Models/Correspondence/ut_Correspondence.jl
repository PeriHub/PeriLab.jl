# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using LinearAlgebra: I, inv, norm
using StaticArrays: SMatrix, SVector

# Rigid rotation F = R of a prestressed node with a full neighborhood: the bond forces
# have to rotate with the body (no internal force change in the body frame) and the
# family of the node carries no net moment, sum_j y_j × f_j = 0 (P Fᵀ symmetric).
const CorrespondenceUT = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence
const ZeroEnergyUT = PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control

moment(y::SVector{2}, f::SVector{2}) = y[1] * f[2] - y[2] * f[1]
function moment(y::SVector{3}, f::SVector{3})
    SVector(y[2] * f[3] - y[3] * f[2],
            y[3] * f[1] - y[1] * f[3],
            y[1] * f[2] - y[2] * f[1])
end

for dof in (2, 3)
    @testset "ut_rigid_rotation_prestressed_$(dof)d" begin
        offsets = [SVector{dof,Float64}(Tuple(c))
                   for c in CartesianIndices(ntuple(_ -> -2:2, dof))
                   if 0 < norm(Tuple(c)) <= 2.1]
        nb = length(offsets)
        PeriLab.Data_Manager.initialize_data()
        PeriLab.Data_Manager.set_num_controller(1)
        PeriLab.Data_Manager.set_dof(dof)
        nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                    Int64)
        nn[1] = nb
        D = PeriLab.Data_Manager
        F = D.create_constant_node_tensor_field("Deformation Gradient", Float64, dof)
        Kinv = D.create_constant_node_tensor_field("Inverse Shape Tensor", Float64, dof)
        stress = D.create_constant_node_tensor_field("Cauchy Stress", Float64, dof)
        zstiff = D.create_constant_node_tensor_field("Zero Energy Stiffness", Float64, dof)
        X = D.create_constant_bond_vector_state("Bond Geometry", Float64, dof)
        x = D.create_constant_bond_vector_state("Deformed Bond Geometry", Float64, dof)
        omega = D.create_constant_bond_scalar_state("Influence Function", Float64;
                                                    default_value = 1.0)
        damage = D.create_constant_bond_scalar_state("Bond Damage", Float64;
                                                     default_value = 1.0)
        force = D.create_constant_bond_vector_state("Bond Forces", Float64, dof)
        for j in 1:nb
            X[1][j] = offsets[j]
        end
        Kinv[1, :, :] = inv(sum(xi * xi' for xi in offsets))
        S = dof == 2 ? [300.0 20.0; 20.0 -50.0] :
            [300.0 20.0 -10.0; 20.0 -50.0 15.0; -10.0 15.0 80.0]
        stress[1, :, :] = S
        zstiff[1, :, :] = S   # any stiffness: the control residual F X - x is zero
        R = SMatrix{dof,dof}(PeriLab.Geometry.rotation_tensor(dof == 2 ? [0.03 * 180 / pi] :
                                                              [1.0, 1.5, 0.03 * 180 / pi],
                                                              dof))

        function bond_forces(Fmat)
            F[1, :, :] = Fmat
            for j in 1:nb
                x[1][j] = Fmat * offsets[j]
            end
            CorrespondenceUT.calculate_bond_force!([1], dof, F, X, damage, omega, Kinv,
                                                   stress, force)
            ZeroEnergyUT.zero_energy_mode_forces!([1], zstiff, F, X, x, omega, force)
            return copy(force[1])
        end
        f_ref = bond_forces(Matrix{Float64}(I, dof, dof))
        f_rot = bond_forces(Matrix(R))
        scale = maximum(norm, f_ref)
        @test maximum(norm(f_rot[j] - R * f_ref[j]) for j in 1:nb) < 1e-12 * scale
        @test norm(sum(moment(R * offsets[j], f_rot[j]) for j in 1:nb)) < 1e-12 * scale
        @test norm(sum(f_rot) - R * sum(f_ref)) < 1e-12 * scale
    end
end

# f_ij = d_ij ω_ij F S K⁻¹ ξ_ij for a general (unsymmetric) F, per-bond damage and
# influence, nodes with different neighbor counts and a node with singular F (no bonds).
@testset "ut_calculate_bond_force_values" begin
    # hand computed 2D value: F S = [334 12; 33 -44], F S K⁻¹ = [168.2 36.4; 12.1 -7.7],
    # times ξ = (1, 2) and d ω = 0.4 gives (96.4, -1.32)
    for dof in (2, 3)
        D = PeriLab.Data_Manager
        D.initialize_data()
        D.set_num_controller(3)
        D.set_dof(dof)
        D.create_constant_node_scalar_field("Number of Neighbors", Int64) .= [2, 3, 1]
        F = D.create_constant_node_tensor_field("Deformation Gradient", Float64, dof)
        Kinv = D.create_constant_node_tensor_field("Inverse Shape Tensor", Float64, dof)
        stress = D.create_constant_node_tensor_field("Cauchy Stress", Float64, dof)
        X = D.create_constant_bond_vector_state("Bond Geometry", Float64, dof)
        omega = D.create_constant_bond_scalar_state("Influence Function", Float64)
        damage = D.create_constant_bond_scalar_state("Bond Damage", Float64)
        force = D.create_constant_bond_vector_state("Bond Forces", Float64, dof)
        if dof == 2
            Fs = ([1.1 0.2; 0.05 0.9], [0.95 -0.1; 0.3 1.2])
            Ss = ([300.0 20.0; 20.0 -50.0], [-80.0 35.0; 35.0 140.0])
            Ks = ([0.5 0.1; 0.1 0.25], [0.7 -0.2; -0.2 0.4])
        else
            Fs = ([1.1 0.2 -0.1; 0.05 0.9 0.15; 0.0 -0.2 1.05],
                  [0.95 -0.1 0.2; 0.3 1.2 0.0; -0.15 0.1 0.85])
            Ss = ([300.0 20.0 -10.0; 20.0 -50.0 15.0; -10.0 15.0 80.0],
                  [-80.0 35.0 5.0; 35.0 140.0 -25.0; 5.0 -25.0 60.0])
            Ks = ([0.5 0.1 0.0; 0.1 0.25 0.05; 0.0 0.05 0.3],
                  [0.7 -0.2 0.1; -0.2 0.4 0.0; 0.1 0.0 0.6])
        end
        for i in 1:2
            F[i, :, :] = Fs[i]
            stress[i, :, :] = Ss[i]
            Kinv[i, :, :] = Ks[i]
        end
        F[3, :, :] .= 0.0                    # node without bonds
        stress[3, :, :] = Ss[1]
        Kinv[3, :, :] = Ks[1]
        xis = [SVector{dof,Float64}(ntuple(k -> (j + k)^2 / 4 - k, dof)) for j in 1:6]
        b = 0
        for i in 1:3, j in eachindex(X[i])
            b += 1
            X[i][j] = xis[b]
            omega[i][j] = 0.5 + 0.1 * b
            damage[i][j] = isodd(b) ? 0.8 : 1.0
        end
        X[1][1] = SVector{dof,Float64}(ntuple(k -> Float64(k), dof))   # ξ = (1, 2[, 3])
        omega[1][1], damage[1][1] = 0.5, 0.8

        CorrespondenceUT.calculate_bond_force!([1, 2, 3], dof, F, X, damage, omega, Kinv,
                                               stress, force)
        dof == 2 && @test force[1][1] ≈ SVector(96.4, -1.32) rtol = 1e-14
        for i in 1:2, j in eachindex(X[i])
            P = Fs[i] * Ss[i]
            @test force[i][j] ≈ damage[i][j] * omega[i][j] * (P * Ks[i] * X[i][j]) rtol = 1e-13
        end
        @test all(iszero, force[3])
    end
end
