# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

@testset "contact_initialize_data" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    params = Dict()
    PeriLab.Solver_Manager.Model_Factory.Contact.Penalty_Model.init_contact_model(params)
    @test params["Contact Stiffness"] == 1e8
    params["Friction Coefficient"] = -3
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Contact.Penalty_Model.init_contact_model(params))
end
