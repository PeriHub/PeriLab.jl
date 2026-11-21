# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const CEA = PeriLab.Solver_Manager.Model_Factory.Damage.Critical_Energy_Aniso
const CEA_DAM = PeriLab.Solver_Manager.Model_Factory.Damage

# Two nodes, one bond pair at the global angle `bond_angle` (degrees), material axes rotated by
# `phi`; the bond is stretched by `stretch` and carries the half bond energy `energy`.
function ut_aniso_setup(raw; bond_angle = 0.0, phi = 0.0, stretch = 0.1, energy = 1.0)
    DM = PeriLab.Data_Manager
    DM.initialize_data()
    DM.set_num_controller(2)
    DM.set_dof(2)
    DM.create_constant_node_scalar_field("Number of Neighbors", Int64) .= 1
    nlist = DM.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] = [2]
    nlist[2] = [1]
    DM.set_inverse_nlist(PeriLab.Helpers.find_inverse_bond_id(nlist))
    DM.create_constant_node_scalar_field("Block_Id", Int64) .= 1
    DM.create_constant_node_scalar_field("Update", Bool)
    DM.create_constant_node_scalar_field("Horizon", Float64) .= 1.0
    DM.create_constant_node_scalar_field("Volume", Float64) .= 1.0
    DM.create_constant_node_scalar_field("Dissipated Energy", Float64)
    DM.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    len = DM.create_constant_bond_scalar_state("Bond Length", Float64)
    _, dlen = DM.create_bond_scalar_state("Deformed Bond Length", Float64)
    for i in 1:2
        len[i] .= 1.0
        dlen[i] .= 1 + stretch
    end
    S = DM.SVector{2,Float64}
    e = S(cosd(bond_angle), sind(bond_angle))
    # w = 0.25 |(f_ij - f_ji) . (y - x)| = 0.25 * 2c * stretch
    c = 2 * energy / stretch
    bg = DM.create_constant_bond_vector_state("Bond Geometry", Float64, 2)
    _, dbg = DM.create_bond_vector_state("Deformed Bond Geometry", Float64, 2)
    forces = DM.create_constant_bond_vector_state("Bond Forces", Float64, 2)
    for (i, sign) in ((1, 1.0), (2, -1.0))
        bg[i][1] = sign * e
        dbg[i][1] = sign * (1 + stretch) * e
        forces[i][1] = sign * c * e
    end
    R = DM.create_constant_node_tensor_field("Rotation Tensor", Float64, 2)
    for iID in 1:2
        R[iID, :, :] = [cosd(phi) -sind(phi); sind(phi) cosd(phi)]
    end
    DM.set_rotation(true)
    d = ut_typed_damage(merge(Dict{String,Any}("Damage Model" => "Critical Energy Anisotropic",
                                               "Critical Value" => 1.0), raw))
    CEA_DAM.init_aniso_crit_values(d.base.anisotropic_damage, 1, 2)
    CEA.init_model([1, 2], d.model, d, 1)
    DM.get_field("Quad Horizon") .= 1.0     # critical bond energy = critical value
    return d
end

function ut_aniso_damage(raw; kwargs...)
    d = ut_aniso_setup(raw; kwargs...)
    CEA.compute_model([1, 2], d.model, d, 1, 0.0, 1.0)
    bd = PeriLab.Data_Manager.get_bond_damage("NP1")
    @test bd[1][1] ≈ bd[2][1]    # both half bonds see the same weights
    return bd[1][1]
end

ut_scale_forces!(s) = foreach(f -> (f .*= s), PeriLab.Data_Manager.get_field("Bond Forces"))

function ut_aniso(criterion; x = 2.0, y = 0.5, extra...)
    Dict{String,Any}("Anisotropic Damage" => Dict{String,
                          Any}("Critical Value X" => x,
                               "Critical Value Y" => y,
                               "Bond Criterion" => criterion,
                               (string(k) => v
                                for (k, v) in extra)...))
end

@testset "anisotropic Axis criterion" begin
    # G_x = 2, G_y = 0.5, w = 1: only the y axis fails, the bond keeps the x share n_x^2
    p = ut_aniso("Axis")
    @test ut_aniso_damage(p; bond_angle = 0.0) == 1.0
    @test ut_aniso_damage(p; bond_angle = 90.0) == 0.0
    @test ut_aniso_damage(p; bond_angle = 45.0) ≈ 0.5
    @test ut_aniso_damage(p; bond_angle = -45.0) ≈ 0.5
    @test ut_aniso_damage(p; bond_angle = 30.0) ≈ 0.75
    # rotated material, phi = 30: the material x axis is the global 30 degree direction
    @test ut_aniso_damage(p; bond_angle = 30.0, phi = 30.0) ≈ 1.0
    @test ut_aniso_damage(p; bond_angle = -30.0, phi = 30.0) ≈ 0.25
    @test ut_aniso_damage(p; bond_angle = 120.0, phi = 30.0) ≈ 0.0 atol=1e-12
    # w above both critical values
    @test ut_aniso_damage(p; bond_angle = 30.0, energy = 3.0) == 0.0
end

@testset "anisotropic Axis irreversibility and dissipated energy" begin
    DM = PeriLab.Data_Manager
    d = ut_aniso_setup(ut_aniso("Axis"); bond_angle = 45.0)
    CEA.compute_model([1, 2], d.model, d, 1, 0.0, 1.0)
    @test DM.get_bond_damage("NP1")[1][1] ≈ 0.5
    @test DM.get_field("Anisotropic Axis Failure")[1][1] == [false, true]
    @test DM.get_field("Update") == [true, true]
    # w (1) * lost share (0.5) * V_i V_j
    @test DM.get_field("Dissipated Energy") ≈ [0.5, 0.5]
    # unloading: the failed axis stays failed, nothing changes
    DM.get_field("Update") .= false
    ut_scale_forces!(0.1)
    CEA.compute_model([1, 2], d.model, d, 1, 0.0, 1.0)
    @test DM.get_bond_damage("NP1")[1][1] ≈ 0.5
    @test DM.get_field("Anisotropic Axis Failure")[1][1] == [false, true]
    @test DM.get_field("Update") == [false, false]
    @test DM.get_field("Dissipated Energy") ≈ [0.5, 0.5]
    # w = 3 > G_x: the x axis fails too; the remaining half is released
    ut_scale_forces!(30.0)
    CEA.compute_model([1, 2], d.model, d, 1, 0.0, 1.0)
    @test DM.get_bond_damage("NP1")[1][1] == 0.0
    @test DM.get_field("Anisotropic Axis Failure")[1][1] == [true, true]
    @test DM.get_field("Dissipated Energy") ≈ [0.5 + 3 * 1.0, 0.5 + 3 * 1.0]
    @test !DM.has_key("Anisotropic Force Filter")
end

@testset "anisotropic Axis Force criterion (dominance rule)" begin
    p = ut_aniso("Axis Force")
    @test ut_aniso_damage(p; bond_angle = 30.0) ≈ 0.75      # w_y = 0.25 < 0.5
    @test ut_aniso_damage(p; bond_angle = 45.0) == 0.0      # w_y = 0.5: dominant
    @test ut_aniso_damage(p; bond_angle = -45.0) == 0.0
    @test ut_aniso_damage(p; bond_angle = 30.0, phi = 30.0) ≈ 1.0
    @test ut_aniso_damage(p; bond_angle = -30.0, phi = 30.0) == 0.0  # w_y = 0.75
    ut_aniso_setup(p)
    @test PeriLab.Data_Manager.has_key("Anisotropic Force Filter")
    ut_aniso_setup(ut_aniso("Axis Force"; var"Force Filter" = false))
    @test !PeriLab.Data_Manager.has_key("Anisotropic Force Filter")
end

@testset "anisotropic Directional criterion" begin
    # 45 degrees: G = (0.5/2 + 0.5/0.5)^-1 = 0.8
    p = ut_aniso("Directional")
    @test ut_aniso_damage(p; bond_angle = 45.0, energy = 0.9) == 0.0
    @test ut_aniso_damage(p; bond_angle = -45.0, energy = 0.9) == 0.0
    @test ut_aniso_damage(p; bond_angle = 45.0, energy = 0.7) == 1.0
    @test ut_aniso_damage(p; bond_angle = 0.0, energy = 1.9) == 1.0
    @test ut_aniso_damage(p; bond_angle = 30.0, phi = 30.0, energy = 1.9) == 1.0
    @test ut_aniso_damage(p; bond_angle = -60.0, phi = 30.0, energy = 1.9) == 0.0
    ut_aniso_damage(p; bond_angle = 45.0, energy = 0.9)
    @test PeriLab.Data_Manager.get_field("Dissipated Energy") ≈ [0.9, 0.9]
end

@testset "anisotropic Axis Stretch criterion" begin
    # critical stretches s_x = 0.2, s_y = 0.05, bond stretch 0.1
    p = ut_aniso("Axis Stretch"; x = 0.2, y = 0.05)
    @test ut_aniso_damage(p; bond_angle = 0.0) == 1.0
    @test ut_aniso_damage(p; bond_angle = 30.0) ≈ 0.75
    @test ut_aniso_damage(p; bond_angle = 45.0) == 0.0
    @test ut_aniso_damage(p; bond_angle = -45.0) == 0.0
    @test ut_aniso_damage(p; bond_angle = 30.0, phi = 30.0) ≈ 1.0
    @test ut_aniso_damage(p; bond_angle = 0.0, phi = 30.0) ≈ 0.75
    @test ut_aniso_damage(p; bond_angle = 0.0, stretch = 0.3) == 0.0
end

@testset "anisotropic Directional Stretch criterion" begin
    # 45 degrees, n = 1: s_c^2 = 0.05^2 + (0.2^2 - 0.05^2) sqrt(0.5), s_c = 0.1704
    p = ut_aniso("Directional Stretch"; x = 0.2, y = 0.05)
    @test ut_aniso_damage(p; bond_angle = 45.0, stretch = 0.18) == 0.0
    @test ut_aniso_damage(p; bond_angle = -45.0, stretch = 0.18) == 0.0
    @test ut_aniso_damage(p; bond_angle = 45.0, stretch = 0.16) == 1.0
    @test ut_aniso_damage(p; bond_angle = 30.0, phi = 30.0, stretch = 0.19) == 1.0
    @test ut_aniso_damage(p; bond_angle = -30.0, phi = 30.0, stretch = 0.19) == 0.0
end

@testset "anisotropic force filter" begin
    DM = PeriLab.Data_Manager
    COR = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence
    ut_aniso_setup(ut_aniso("Axis Force"); phi = 30.0)
    failure = DM.get_field("Anisotropic Axis Failure")
    failure[1][1] = DM.SVector(false, true)
    forces = DM.get_field("Bond Forces")
    forces[1][1] = DM.SVector(1.0, 2.0)
    forces[2][1] = DM.SVector(1.0, 2.0)
    COR.filter_failed_axes!([1, 2], forces, DM.get_bond_damage("NP1"))
    a = DM.SVector(cosd(30), sind(30))      # the intact material x axis
    @test forces[1][1] ≈ (a' * DM.SVector(1.0, 2.0)) * a
    @test forces[2][1] == DM.SVector(1.0, 2.0)
end

@testset "get_quad_horizon anisotropic" begin
    @test CEA.get_quad_horizon(2.0, 2, 3.0) == 3 / (4 * 2.0^3 * 3.0)
    @test CEA.get_quad_horizon(2.0, 3, 3.0) == 2 / (pi * 2.0^4)
end
