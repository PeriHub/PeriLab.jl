# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Models/Contact/Contact_search.jl")
# include("../../../../src/Core/Data_manager.jl")
using Test
using .Contact_search: init_contact_search,
                       find_potential_contact_pairs, create_potential_contact_dict

@testset "ut_find_potential_contact_pairs and ut_create_potential_contact_dict" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(2)
    points = [0.0 0.0;
              2.32241 2.22618;
              0.17273 4.04840;
              2.74944 3.88023;
              6.44488 3.57078]
    test_data_manager.set_free_surface_nodes(1, [1, 4, 5])
    test_data_manager.set_free_surface_nodes(3, [2, 3])
    test_data_manager.set_all_positions(points)
    contact_params = Dict("Master Block ID" => 3, "Slave Block ID" => 1,
                          "Contact Radius" => 0.01)
    nearpoints = find_potential_contact_pairs(test_data_manager, contact_params)
    @test nearpoints == [[], []]
    @test create_potential_contact_dict(nearpoints, test_data_manager, contact_params) ==
          Dict{Int64,Vector{Int64}}()

    contact_params = Dict("Master Block ID" => 3, "Slave Block ID" => 1,
                          "Contact Radius" => 2)
    nearpoints = find_potential_contact_pairs(test_data_manager, contact_params)
    @test nearpoints == [[2], []]
    @test create_potential_contact_dict(nearpoints, test_data_manager, contact_params) ==
          Dict{Int64,Vector{Int64}}(2 => [4])

    contact_params = Dict("Master Block ID" => 3, "Slave Block ID" => 1,
                          "Contact Radius" => 4)
    nearpoints = find_potential_contact_pairs(test_data_manager, contact_params)
    @test nearpoints == [[1, 2], [2]]
    @test create_potential_contact_dict(nearpoints, test_data_manager, contact_params) ==
          Dict{Int64,Vector{Int64}}(2 => [1, 4], 3 => [4])
    contact_params = Dict("Master Block ID" => 3, "Slave Block ID" => 1,
                          "Contact Radius" => 8)
    nearpoints = find_potential_contact_pairs(test_data_manager, contact_params)
    @test nearpoints == [[1, 2, 3], [1, 2, 3]]
    @test create_potential_contact_dict(nearpoints, test_data_manager, contact_params) ==
          Dict{Int64,Vector{Int64}}(2 => [1, 4, 5], 3 => [1, 4, 5])
end
