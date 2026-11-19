# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
# include("../../../../../../src/PeriLab.jl")
# using .PeriLab

@testset "get_name&fe_support" begin
    @test parentmodule(PeriLab.ParameterSpec.lookup_model(:material,
                                                          "Correspondence Plastic")) ===
          PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic
    @test !(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.fe_support())
end

@testset "ut_init_model" begin
    nodes = 2
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(nodes)
    PeriLab.Data_Manager.set_dof(3)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn .= 2
    PLASTIC = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                  "Shear Modulus" => 10.5, "Yield Stress" => 3.4))
    PLASTIC.init_model(Vector{Int64}(1:nodes), m.model.parts[2], m, 1)
    @test PeriLab.Data_Manager.has_key("von Mises Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic StrainN")
    ba = typed_block_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                                   "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                   "Shear Modulus" => 10.5, "Yield Stress" => 3.4,
                                   "Bond Associated" => true))
    PLASTIC.init_model(Vector{Int64}(1:nodes), ba.model.parts[2], ba, 1)
    @test PeriLab.Data_Manager.has_key("von Mises Bond Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic Bond StrainN")
end

# one material point under a uniaxial strain path through Correspondence Elastic and
# Correspondence Plastic, as the composite model runs them
function plastic_strain_path(raw; steps = 60, dstrain = 5.0e-4)
    DM = PeriLab.Data_Manager
    DM.initialize_data()
    DM.set_num_controller(1)
    DM.set_dof(2)
    DM.create_constant_node_vector_field("Coordinates", Float64, 2)
    CORR = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence
    m = typed_block_material(raw)
    elastic, plastic = m.model.parts
    CORR.Correspondence_Elastic.init_model([1], elastic, m, 1)
    CORR.Correspondence_Plastic.init_model([1], plastic, m, 1)
    stress_N = zeros(1, 2, 2)
    stress_NP1 = zeros(1, 2, 2)
    inc = zeros(1, 2, 2)
    # uniaxial stress while elastic: the lateral strain keeps sigma_yy = 0
    C = DM.get_field("Material Gradient")
    inc[1, 1, 1] = dstrain
    inc[1, 2, 2] = -C[1, 2, 1] / C[1, 2, 2] * dstrain
    eps_p_N = DM.get_field("Plastic Strain", "N")
    eps_p = DM.get_field("Plastic Strain", "NP1")
    vm = DM.get_field("von Mises Yield Stress", "NP1")
    sxx, ep, vms = Float64[], Float64[], Float64[]
    for _ in 1:steps
        CORR.Correspondence_Elastic.compute_stresses([1], 2, elastic, m, 0.0, 1.0, inc,
                                                     stress_N, stress_NP1)
        CORR.Correspondence_Plastic.compute_stresses([1], 2, plastic, m, 0.0, 1.0, inc,
                                                     stress_N, stress_NP1)
        push!(sxx, stress_NP1[1, 1, 1])
        push!(ep, eps_p[1])
        push!(vms, vm[1])
        stress_N .= stress_NP1
        eps_p_N .= eps_p
    end
    return (; sxx, ep, vms, m)
end

@testset "ut_plastic_anisotropic_elasticity" begin
    E, nu, sy = 5836.04, 0.33, 75.2
    plastic = "Correspondence Elastic + Correspondence Plastic"
    iso = plastic_strain_path(Dict("Material Model" => plastic,
                                   "Symmetry" => "isotropic plane stress",
                                   "Young's Modulus" => E, "Poisson's Ratio" => nu,
                                   "Yield Stress" => sy))
    @test maximum(iso.ep) > 0   # the path goes past yield
    # the return mapping keeps the hydrostatic stress: no drop at yield
    @test all(diff(iso.sxx) .> -1e-9)
    # transverse isotropic, isotropic in disguise: same curve and plastic strain
    disguise = plastic_strain_path(Dict("Material Model" => plastic,
                                        "Symmetry" => "transverse isotropic plane stress",
                                        "Young's Modulus X" => E,
                                        "Young's Modulus Y" => E,
                                        "Shear Modulus XY" => E / (2 * (1 + nu)),
                                        "Poisson's Ratio XY" => nu,
                                        "Yield Stress" => sy))
    @test disguise.m.moduli === nothing
    @test disguise.sxx ≈ iso.sxx
    @test disguise.ep ≈ iso.ep

    # truly transverse isotropic (ENF adhesive), with the extra Shear Modulus key
    raw = Dict("Material Model" => plastic,
               "Symmetry" => "transverse isotropic plane stress",
               "Young's Modulus X" => E, "Young's Modulus Y" => 13844.0,
               "Shear Modulus XY" => 2193.94, "Poisson's Ratio XY" => nu,
               "Yield Stress" => sy)
    ti = plastic_strain_path(raw)
    @test PeriLab.Data_Manager.get_field("Plastic Shear Modulus")[1] == 2193.94
    yielded = ti.ep .> 0
    @test any(yielded) && !all(yielded)
    @test all(ti.vms[.!yielded] .< sy)
    @test all(ti.vms[yielded] .≈ sy)
    @test all(diff(ti.ep) .>= 0) && ti.ep[end] > ti.ep[findfirst(yielded)]
    plastic_strain_path(merge(raw, Dict("Shear Modulus" => 2000.0)))
    @test PeriLab.Data_Manager.get_field("Plastic Shear Modulus")[1] == 2000.0

    # full stiffness matrix: the in-plane shear entry C66
    C = Dict("C$i$j" =>
                 (i == j ? (i <= 3 ? 10000.0 : 3000.0) : (i <= 3 && j <= 3 ? 3000.0 : 0.0))
             for i in 1:6 for j in i:6)
    plastic_strain_path(merge(C,
                              Dict("Material Model" => plastic,
                                   "Symmetry" => "anisotropic plane strain",
                                   "Yield Stress" => sy)); steps = 1)
    @test PeriLab.Data_Manager.get_field("Plastic Shear Modulus")[1] == 3000.0
end
