# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

@testset "ut_check_valid_contact_model" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(2)
    test_data_manager.set_num_controller(4)
    block_id = test_data_manager.create_constant_node_field("Block_Id", Int64, 1)
    block_id .= 1
    # contact_params = Dict("cm" => Dict("Contact Groups" => Dict()))
    # @test !PeriLab.Solver_Manager.check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Contact Groups" => Dict("cg" => Dict("Master Block ID" => 1))))
    @test !PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Contact Groups" => Dict("cg" => Dict("Master Block ID" => 1,
                                                                             "Slave Block ID" => 1))))
    @test !PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Contact Groups" => Dict("cg" => Dict("Master Block ID" => 1,
                                                                             "Slave Block ID" => 2))))
    @test !PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Contact Groups" => Dict("cg" => Dict("Master Block ID" => 1,
                                                                             "Slave Block ID" => 2))),
                          "cm2" => Dict("Contact Groups" => Dict("cg" => Dict("Master Block ID" => 2,
                                                                              "Slave Block ID" => 1))))
    @test !PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact_params, block_id)
    block_id[2] = 2
    @test !PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Contact Groups" => Dict("cg" => Dict("Master Block ID" => 1,
                                                                             "Slave Block ID" => 2,
                                                                             "Search Radius" => 0.0))))
    @test !PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact_params, block_id)
    contact_params = Dict("cm" => Dict("Contact Groups" => Dict("cg" => Dict("Master Block ID" => 1,
                                                                             "Slave Block ID" => 2,
                                                                             "Search Radius" => -20.0))))
    @test !PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact_params, block_id)
end

@testset "ut_get_all_contact_blocks" begin
    params = Dict("cm" => Dict("Contact Groups" => Dict("a" => Dict("Master Block ID" => 1,
                                                                    "Slave Block ID" => 2),
                                                        "ba" => Dict("Master Block ID" => 2,
                                                                     "Slave Block ID" => 1),
                                                        "c" => Dict("Master Block ID" => 3,
                                                                    "Slave Block ID" => 5),
                                                        "q" => Dict("Master Block ID" => 8,
                                                                    "Slave Block ID" => 5))))
    @test PeriLab.Solver_Manager.Model_Factory.Contact.get_all_contact_blocks(params) == [1, 2, 3, 5, 8]
end
@testset "ut_get_double_surfs" begin
    @warn "TBD implentation of double surfs test"
    #get_double_surfs(normals_i, offsets_i, normals_j, offsets_j)
end
