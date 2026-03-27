# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
using TimerOutputs
#include("../../../../../../src/PeriLab.jl")
# #using .PeriLab

# const to = TimerOutput()

@testset "material_name" begin
    @test PeriLab.Solver_Manager.Model_Factory.Material.OneD_Bond_Based_Elastic.material_name() ==
          "1D Bond-based Elastic"
    @test !(PeriLab.Solver_Manager.Model_Factory.Material.OneD_Bond_Based_Elastic.fe_support())
end
@testset "compute_model" begin
    nodes = 2

    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(nodes)
    dof = 3
    PeriLab.Data_Manager.set_dof(dof)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 2
    nn[2] = 3
    h = PeriLab.Data_Manager.create_constant_node_scalar_field("Horizon", Float64)

    h[1:nodes] = 1:nodes
    bf = PeriLab.Data_Manager.create_constant_bond_vector_state("Bond Forces", Float64, dof)

    bdN,
    bdNP1 = PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64;
                                                          default_value = 1)
    dbN,
    dbNP1 = PeriLab.Data_Manager.create_bond_vector_state("Deformed Bond Geometry", Float64,
                                                          dof; default_value = 1)
    dbdN,
    dbdNP1 = PeriLab.Data_Manager.create_bond_scalar_state("Deformed Bond Length", Float64)
    bg = PeriLab.Data_Manager.create_bond_vector_state("Bond Geometry", Float64, dof)
    bd = PeriLab.Data_Manager.create_constant_bond_scalar_state("Bond Length", Float64;
                                                                default_value = 1)
    for iID in 1:nodes
        dbdNP1[iID] .= 1 + (-1)^iID * 0.1
    end
    PeriLab.Solver_Manager.Model_Factory.Material.OneD_Bond_Based_Elastic.init_model(Vector{Int64}(1:nodes),
                                                                                     Dict("Bulk Modulus" => 1.0,
                                                                                          "Young's Modulus" => 1.0))
    # PeriLab.Solver_Manager.Model_Factory.Material.OneD_Bond_Based_Elastic.compute_model(Vector{Int64}(1:nodes),
    #                                                                                     Dict("Bulk Modulus" => 1.0,
    #                                                                                          "Young's Modulus" => 1.0,
    #                                                                                          "Id1" => 1,
    #                                                                                          "Id2" => 2),
    #                                                                                     1,
    #                                                                                     0.0,
    #                                                                                     0.0)
end
