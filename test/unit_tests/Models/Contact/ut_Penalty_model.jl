# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../src/Models/Contact/Penalty_model.jl")

using Test

@testset "contact_initialize_data" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    params = Dict()
    Penalty_model.init_contact_model(test_data_manager, params)
    @test params["Contact Stiffness"] == 1e8
    params["Friction Coefficient"] = -3
    @test isnothing(Penalty_model.init_contact_model(test_data_manager, params))
end
