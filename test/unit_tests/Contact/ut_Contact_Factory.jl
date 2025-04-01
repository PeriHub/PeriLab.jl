# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Models/Contact/Contact_Factory.jl")
# include("../../../../src/Core/Data_manager.jl")
using Test
using .Contact_Factory: init_contact_model

@testset "ut_init_contact_exceptions" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(2)
    test_data_manager.set_num_controller(4)
    block_id = test_data_manager.create_field("Block_Id", Int64, 1)
    block_id .= 1
    contact_params = Dict()
    @test isnothing(init_contact_model(datamanager, contact_params))
    contact_params = Dict("Master" => 1)
    @test isnothing(init_contact_model(datamanager, contact_params))
    contact_params = Dict("Master" => 1, "Slave" => 1)
    @test isnothing(init_contact_model(datamanager, contact_params))
    contact_params = Dict("Master" => 1, "Slave" => 2)
    @test isnothing(init_contact_model(datamanager, contact_params))
    block_id[2] = 2
    @test isnothing(init_contact_model(datamanager, contact_params))
end
