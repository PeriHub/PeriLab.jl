# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

@testset "init_fields" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(3)
    PeriLab.Data_Manager.set_num_controller(4)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist",
                                                                   Int64)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1]
    nlist[4] = [1, 3]
    PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64)
    @test_logs (:error,
                "'Activation_Time' is missing. Please define an 'Activation_Time' for each point in the mesh file.") PeriLab.Solver_Manager.Model_Factory.Additive.init_fields()
    PeriLab.Data_Manager.create_constant_node_scalar_field("Activation_Time", Float64)
    PeriLab.Solver_Manager.Model_Factory.Additive.init_fields()
    field_keys = PeriLab.Data_Manager.get_all_field_keys()
    @test "Active" in field_keys
    active = PeriLab.Data_Manager.get_field("Active")
    for iID in 1:4
        @test active[iID] == false
    end
end

@testset "init_additive" begin
    PeriLab.Data_Manager.data["properties"][23] = Dict("Additive Model" => Dict("Additive Model" => "does not exist"))
    @test_logs (:error, "No additive model of name does not exist exists.") PeriLab.Solver_Manager.Model_Factory.Additive.init_model(Vector{Int64}([
                                                                                                                                                       1,
                                                                                                                                                       2,
                                                                                                                                                       3
                                                                                                                                                   ]),
                                                                                                                                     23)
end
