# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const CMID = PeriLab.InputDeck

function ut_category_deck(models; block = Dict{String,Any}())
    return Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("block_1" => merge(Dict{String,Any}("Block ID" => 1,
                                                                                             "Density" => 1.0,
                                                                                             "Horizon" => 1.0),
                                                                            block)),
                            "Models" => Dict{String,Any}(models),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
end
ut_category_messages(ctx) = Dict(e.path => e.message for e in ctx.errors)

@testset "additive input errors" begin
    _, ctx = CMID.read_input(ut_category_deck(Dict("Additive Models" => Dict("Add" => Dict("Additive Model" => "Simple",
                                                                                             "Print Temperature" => 100.0)))))
    msg = ut_category_messages(ctx)["Models.\"Additive Models\".Add.\"Additive Model\""]
    @test startswith(msg, "model \"Simple\"")
    input, ctx = CMID.read_input(ut_category_deck(Dict{String,Any}()))
    @test isempty(input.additives)
end
