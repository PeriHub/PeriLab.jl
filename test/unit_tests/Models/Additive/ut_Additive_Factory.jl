# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../src/Models/Additive/Additive_Factory.jl")
# include("../../../../src/Core/data_manager.jl")


using Test
using .Additive



@testset "init_fields" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(4)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1]
    nlist[4] = [1, 3]
    test_data_manager.create_bond_field("Bond Damage", Float64, 1)
    test_data_manager = Additive.init_fields(test_data_manager)
    field_keys = test_data_manager.get_all_field_keys()
    @test "Active List" in field_keys
    active = test_data_manager.get_field("Active List")
    for iID = 1:4
        @test active[iID] == false
    end
end

@testset "init_additive" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.data["properties"][23] =
        Dict("Additive Model" => Dict("Additive Model" => "does not exist"))
    @test isnothing(Additive.init_model(test_data_manager, Vector{Int64}([1, 2, 3]), 23))
end
