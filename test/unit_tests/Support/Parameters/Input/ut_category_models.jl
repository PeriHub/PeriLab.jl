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

@testset "degradation models are parsed into typed structs" begin
    input, ctx = CMID.read_input(ut_category_deck(Dict("Degradation Models" => Dict("Deg" => Dict("Degradation Model" => "Thermal Decomposition",
                                                                                                    "Decomposition Temperature" => 50)))))
    @test isempty(ctx.errors)
    @test input.degradations["Deg"].decomposition_temperature == 50.0
    _, ctx = CMID.read_input(ut_category_deck(Dict("Degradation Models" => Dict("Deg" => Dict("Degradation Model" => "Thermal Decomposition")))))
    @test ut_category_messages(ctx)["Models.\"Degradation Models\".Deg.\"Decomposition Temperature\""] ==
          "missing (required by Thermal Decomposition)"
end

function ut_thermal(entry)
    return CMID.read_input(ut_category_deck(Dict("Thermal Models" => Dict("Th" => Dict{String,Any}(entry)))))
end
const TH_PATH = "Models.\"Thermal Models\".Th"

@testset "thermal models are parsed into typed structs" begin
    input, ctx = ut_thermal(Dict("Thermal Model" => "Thermal Flow + Heat Transfer",
                                 "Thermal Conductivity" => 2.0, "Type" => "Correspondence",
                                 "Heat Transfer Coefficient" => 1.0,
                                 "Environmental Temperature" => "20+t"))
    @test isempty(ctx.errors)
    th = input.thermals["Th"]
    @test th.base.thermal_conductivity == 2.0
    flow, transfer = th.model.parts
    @test string(flow.type) == "Correspondence"
    @test transfer.environmental_temperature == "20+t"
    @test transfer.allow_surface_change
    input, _ = ut_thermal(Dict("Thermal Model" => "Thermal Flow+Heat Transfer",
                               "Heat Transfer Coefficient" => 1.0,
                               "Environmental Temperature" => 30))
    flow, transfer = input.thermals["Th"].model.parts
    @test string(flow.type) == "BondBased"                     # default
    @test transfer.environmental_temperature === 30.0
    @test input.thermals["Th"].base.thermal_conductivity === nothing
    input, _ = ut_thermal(Dict("Thermal Model" => "Thermal Expansion",
                               "Thermal Expansion Coefficient" => [1.0, 2.0]))
    @test input.thermals["Th"].model.thermal_expansion_coefficient == [1.0, 2.0]
    @test input.thermals["Th"].model.reference_temperature === nothing
    input, _ = ut_thermal(Dict("Thermal Model" => "HETVAL", "File" => "h.so",
                               "Property_1" => 2.0))
    @test input.thermals["Th"].model.hetval_name == "HETVAL"
    @test input.thermals["Th"].model.number_of_state_variables == 1
end

@testset "thermal input errors" begin
    _, ctx = ut_thermal(Dict("Thermal Model" => "Thermal Flow", "Type" => "Bond"))
    @test haskey(ut_category_messages(ctx), "$TH_PATH.Type")
    _, ctx = ut_thermal(Dict("Thermal Model" => "Thermal Flow", "Thermal Conductivity" => 1.0,
                             "Print Bed Temperature" => 300.0))
    @test ut_category_messages(ctx)["$TH_PATH.\"Thermal Conductivity Print Bed\""] ==
          "required when Print Bed Temperature is given"
    _, ctx = ut_thermal(Dict("Thermal Model" => "HETVAL", "File" => "h.so",
                             "HETVAL Material Name" => repeat("x", 81)))
    @test ut_category_messages(ctx)["$TH_PATH.\"HETVAL Material Name\""] ==
          "at most 80 characters (Fortran)"
    _, ctx = ut_thermal(Dict("Thermal Model" => "Heat Transfer",
                             "Heat Transfer Coefficient" => 1.0,
                             "Environmental Temperature" => 30, "Type" => "Bond based"))
    @test startswith(ut_category_messages(ctx)["$TH_PATH.Type"], "unknown key")
end

@testset "pre-calculation switches" begin
    global_switches = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false,
                           "Bond Associated Deformation Gradient" => false)
    input, ctx = CMID.read_input(ut_category_deck(Dict("Pre Calculation Global" => global_switches,
                                                       "Pre Calculation Models" => Dict("Pre" => Dict("Deformation Gradient" => true)))))
    @test isempty(ctx.errors)
    @test input.pre_calculation_global ==
          Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false)
    @test input.pre_calculations["Pre"] == Dict("Deformation Gradient" => true)
    _, ctx = CMID.read_input(ut_category_deck(Dict("Pre Calculation Global" => Dict("Bond Associated Deformation Gradient" => true))))
    @test ut_category_messages(ctx)["Models.\"Pre Calculation Global\".\"Bond Associated Deformation Gradient\""] ==
          "no longer supported; use \"Bond Associated Correspondence\""
    _, ctx = CMID.read_input(ut_category_deck(Dict("Pre Calculation Global" => Dict("Shape Tenser" => true))))
    @test ut_category_messages(ctx)["Models.\"Pre Calculation Global\".\"Shape Tenser\""] ==
          "unknown key — did you mean \"Shape Tensor\"?"
    input, _ = CMID.read_input(ut_category_deck(Dict{String,Any}()))
    @test input.pre_calculation_global === nothing
    @test isempty(input.pre_calculations)
end
