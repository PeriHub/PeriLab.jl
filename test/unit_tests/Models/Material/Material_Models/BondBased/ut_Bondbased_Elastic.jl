# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using TimerOutputs
#include("../../../../../../src/PeriLab.jl")
# #using .PeriLab

#const to = TimerOutput()

@testset "material_name" begin
    @test PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.material_name() == "Bond-based Elastic"
    @test !(PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.fe_support())
end
@testset "compute_model" begin
    nodes = 2
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nodes)
    dof = 3
    test_data_manager.set_dof(dof)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    h = test_data_manager.create_constant_node_field("Horizon", Float64, 1)

    h[1:nodes] = 1:nodes
    bf = test_data_manager.create_constant_bond_field("Bond Forces", Float64, dof)

    bdN, bdNP1 = test_data_manager.create_bond_field("Bond Damage", Float64, 1, 1)
    dbN, dbNP1 = test_data_manager.create_bond_field("Deformed Bond Geometry", Float64, dof,
                                                     1)
    dbdN, dbdNP1 = test_data_manager.create_bond_field("Deformed Bond Length", Float64, 1)
    bg = test_data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    bd = test_data_manager.create_constant_bond_field("Bond Length", Float64, 1, 1)
    for iID in 1:nodes
        dbdNP1[iID] .= 1 + (-1)^iID * 0.1
    end
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.init_model(test_data_manager,
                                                     Vector{Int64}(1:nodes),
                                                     Dict("Bulk Modulus" => 1.0,
                                                          "Young's Modulus" => 1.0))
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.compute_model(test_data_manager,
                                                        Vector{Int64}(1:nodes),
                                                        Dict("Bulk Modulus" => 1.0,
                                                             "Young's Modulus" => 1.0),
                                                        1,
                                                        0.0,
                                                        0.0,
                                                        to)

    bf = test_data_manager.get_field("Bond Forces")
    @test isapprox(bf[1],
                   [
                       [-0.21220659078919374, -0.21220659078919374, -0.21220659078919374],
                       [-0.21220659078919374, -0.21220659078919374, -0.21220659078919374]
                   ])
    @test isapprox(bf[2],
                   [
                       [0.010851473392629237, 0.010851473392629237, 0.010851473392629237],
                       [0.010851473392629237, 0.010851473392629237, 0.010851473392629237],
                       [0.010851473392629237, 0.010851473392629237, 0.010851473392629237]
                   ])

    bf[1][1] .= 0
    bf[1][2] .= 0
    bf[2][1] .= 0
    bf[2][2] .= 0
    bf[2][3] .= 0
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.init_model(test_data_manager,
                                                     Vector{Int64}(1:nodes),
                                                     Dict("Bulk Modulus" => 1.0,
                                                          "Young's Modulus" => 1.0,
                                                          "Symmetry" => "here is something"))
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.compute_model(test_data_manager,
                                                        Vector{Int64}(1:nodes),
                                                        Dict("Bulk Modulus" => 1.0,
                                                             "Young's Modulus" => 1.0,
                                                             "Symmetry" => "here is something"),
                                                        1,
                                                        0.0,
                                                        0.0,
                                                        to)

    bf = test_data_manager.get_field("Bond Forces")

    @test isapprox(bf[1],
                   [
                       [-0.21220659078919374, -0.21220659078919374, -0.21220659078919374],
                       [-0.21220659078919374, -0.21220659078919374, -0.21220659078919374]
                   ])
    @test isapprox(bf[2],
                   [
                       [0.010851473392629237, 0.010851473392629237, 0.010851473392629237],
                       [0.010851473392629237, 0.010851473392629237, 0.010851473392629237],
                       [0.010851473392629237, 0.010851473392629237, 0.010851473392629237]
                   ])

    bf[1][1] .= 0
    bf[1][2] .= 0
    bf[2][1] .= 0
    bf[2][2] .= 0
    bf[2][3] .= 0
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.init_model(test_data_manager,
                                                     Vector{Int64}(1:nodes),
                                                     Dict("Bulk Modulus" => 1.0,
                                                          "Young's Modulus" => 1.0,
                                                          "Symmetry" => "plane strain"))

    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.compute_model(test_data_manager,
                                                        Vector{Int64}(1:nodes),
                                                        Dict("Bulk Modulus" => 1.0,
                                                             "Young's Modulus" => 1.0,
                                                             "Symmetry" => "plane strain"),
                                                        1,
                                                        0.0,
                                                        0.0,
                                                        to)

    bf = test_data_manager.get_field("Bond Forces")

    @test isapprox(bf[1],
                   [
                       [-0.16976527263135502, -0.16976527263135502, -0.16976527263135502],
                       [-0.16976527263135502, -0.16976527263135502, -0.16976527263135502]
                   ])
    @test isapprox(bf[2],
                   [
                       [0.01736235742820678, 0.01736235742820678, 0.01736235742820678],
                       [0.01736235742820678, 0.01736235742820678, 0.01736235742820678],
                       [0.01736235742820678, 0.01736235742820678, 0.01736235742820678]
                   ])

    bf[1][1] .= 0
    bf[1][2] .= 0
    bf[2][1] .= 0
    bf[2][2] .= 0
    bf[2][3] .= 0
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.init_model(test_data_manager,
                                                     Vector{Int64}(1:nodes),
                                                     Dict("Bulk Modulus" => 1.0,
                                                          "Young's Modulus" => 1.0,
                                                          "Symmetry" => "plane stress"))
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.compute_model(test_data_manager,
                                                        Vector{Int64}(1:nodes),
                                                        Dict("Bulk Modulus" => 1.0,
                                                             "Young's Modulus" => 1.0,
                                                             "Symmetry" => "plane stress"),
                                                        1,
                                                        0.0,
                                                        0.0,
                                                        to)

    bf = test_data_manager.get_field("Bond Forces")
    @test isapprox(bf[1][1][1], -0.15915494309189532)
    @test isapprox(bf[1][1][2], -0.15915494309189532)
    @test isapprox(bf[1][2][1], -0.15915494309189532)
    @test isapprox(bf[1][2][2], -0.15915494309189532)
    @test isapprox(bf[2][1][1], 0.016277210088943856)
    @test isapprox(bf[2][1][2], 0.016277210088943856)
    @test isapprox(bf[2][2][1], 0.016277210088943856)
    @test isapprox(bf[2][2][2], 0.016277210088943856)
end
