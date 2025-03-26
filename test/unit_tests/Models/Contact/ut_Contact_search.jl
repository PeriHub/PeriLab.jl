# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../src/Models/Contact/Contact_search.jl")
# include("../../../../src/Core/Data_manager.jl")
using Test
using .Corrosion: get_surface_normals, get_surface_connectivity, init_contact
@testset "ut_get_surface_normals" begin
    points_2D = [
        [0.0, 0.0],    # Unten links
        [4, 0.0],      # Unten rechts
        [1, 1],
        [4, 2],        # Oben rechts
        [0.0, 2]       # Oben links
    ]

    @test get_surface_normals(points_2D) == [-0.0 -1.0; -1.0 -0.0; -0.0 1.0; 1.0 -0.0]
    points_3D = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 3.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 0.5],
        [2.0, 0.0, 0.5],
        [2.0, 3.0, 0.5],
        [0.0, 3.0, 0.5]
    ]
    @test get_surface_normals(points_3D) == [-0.0 -0.0 -1.0;
           -1.0 -0.0 -0.0;
           -0.0 -1.0 -0.0;
           -0.0 -0.0 2.0;
           -0.0 1.0 -0.0;
           1.0 -0.0 -0.0]
end

@testset "ut_get_surface_normals" begin
    points_2D = [
        [0.0, 0.0],    # Unten links
        [4, 0.0],      # Unten rechts
        [1, 1],
        [0.5, 0.0],    # Unten links
        [4, 2],        # Oben rechts
        [0.0, 2]       # Oben links
    ]
    test = get_surface_connectivity(points_2D)
    @test haskey(test, 4)
    @test test[4] == 1
end

@testset "ut_init_contact_exceptions" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(2)
    test_data_manager.set_num_controller(4)
    block_id = test_data_manager.create_field("Block_Id")
    block_id .= 1
    contact_params = Dict()
    @test isnothing(init_contact(datamanager, contact_params))
    contact_params = Dict("Master" => 1)
    @test isnothing(init_contact(datamanager, contact_params))
    contact_params = Dict("Master" => 1, "Slave" => 1)
    @test isnothing(init_contact(datamanager, contact_params))
    contact_params = Dict("Master" => 2, "Slave" => 1)
    @test isnothing(init_contact(datamanager, contact_params))
    contact_params = Dict("Master" => 1, "Slave" => 2)
    @test isnothing(init_contact(datamanager, contact_params))
    block_id[2] = 2
    @test isnothing(init_contact(datamanager, contact_params))
end
