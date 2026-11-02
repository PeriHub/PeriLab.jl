# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

ut_cf_group(master, slave) = Dict{String,Any}("Master Block ID" => master,
                                              "Slave Block ID" => slave,
                                              "Search Radius" => 0.01)
ut_cf_model(groups) = Dict{String,Any}("Type" => "Penalty Contact", "Contact Radius" => 0.005,
                                       "Contact Groups" => groups)

@testset "ut_check_valid_contact_model" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_num_controller(4)
    block_id = PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64)
    block_id .= 1
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(1, 2)))))
    @test_logs (:error,
                "Block defintion in slave does not exist.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                               block_id)
    end
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(2, 1)))))
    @test_logs (:error,
                "Block defintion in master does not exist.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                               block_id)
    end
    block_id[2] = 2
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(1, 2))),
                                 "cm2" => ut_cf_model(Dict("cg" => ut_cf_group(2, 1)))))
    @test_logs (:error,
                "Master and Slave should be defined in an inverse way, e.g. Master = 1, Slave = 2 in model 1 and Master = 2, Slave = 1 in model 2.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                               block_id)
    end
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(1, 2)))))
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                                          block_id))
end
@testset "ut_get_double_surfs" begin
    @warn "TBD implentation of double surfs test"
    #get_double_surfs(normals_i, offsets_i, normals_j, offsets_j)
end
