# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

@testset "contact_initialize_data" begin
    PeriLab.Data_Manager.initialize_data()
    contact = typed_contact(Dict("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                         "Contact Radius" => 0.005,
                                                         "Contact Groups" => Dict{String,Any}("g" => Dict{String,Any}("Master Block ID" => 2,
                                                                                                                      "Slave Block ID" => 1,
                                                                                                                      "Search Radius" => 0.01)))))
    model = contact.models["C"]
    PM = PeriLab.Solver_Manager.Model_Factory.Contact.Penalty_Model
    @test model.model isa PM.PenaltyContactParams
    @test isnothing(PM.init_contact_model(model.model, model))
    @test model.model.contact_stiffness == 1e8
end
