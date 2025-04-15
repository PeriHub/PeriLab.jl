# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Models/Contact/Contact_Factory.jl")
# include("../../../../src/Core/Data_manager.jl")
using Test
using .Contact_Factory: check_valid_contact_model

@testset "ut_check_valid_contact_model" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(2)
    test_data_manager.set_num_controller(4)
    block_id = test_data_manager.create_field("Block_Id", Int64, 1)
    block_id .= 1
    contact_params = Dict("cm" => Dict())
    @test !check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Master" => 1))
    @test !check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Master" => 1, "Slave" => 1))
    @test !check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Master" => 1, "Slave" => 2))
    @test !check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Master" => 1, "Slave" => 2),
                          "cm2" => Dict("Master" => 2, "Slave" => 1))
    @test !check_valid_contact_model(contact_params, block_id)
    block_id[2] = 2
    @test !check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Master" => 1, "Slave" => 2, "Search Radius" => 0.0))
    @test !check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Master" => 1, "Slave" => 2,
                                       "Search Radius" => -20.0))
    @test !check_valid_contact_model(contact_params, block_id)
end

@testset "ut_check_valid_contact_model" begin
    params = Dict("a" => Dict("Master" => 1, "Slave" => 2),
                  "ba" => Dict("Master" => 2, "Slave" => 1),
                  "c" => Dict("Master" => 3, "Slave" => 5),
                  "q" => Dict("Master" => 8, "Slave" => 5))
    @test get_all_contact_blocks(params) == [1, 2, 3, 5, 8]
end
