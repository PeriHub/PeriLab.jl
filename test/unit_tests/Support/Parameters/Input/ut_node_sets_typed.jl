# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using DataFrames
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const PH = PeriLab.Parameter_Handling

@testset "read_node_sets equals get_node_sets" begin
    dir = mktempdir()
    write(joinpath(dir, "ns.txt"), "header: global_id\n2\n3\n")
    sets = Dict{String,Any}("id" => 2, "list" => "1 3", "range" => "1:3",
                            "expr" => "x > 0.5", "all" => "All", "file" => "ns.txt")
    raw = Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                           "Node Sets" => sets)
    ctx = PS.ParseContext()
    d = PS.convert_value(ID.DiscretizationParams, raw, "D", ctx)
    @test isempty(ctx.errors)
    mesh = DataFrame(x = [0.0, 1.0, 2.0], y = [0.0, 0.0, 0.0])
    typed = PH.read_node_sets(d, dir, mesh)
    reference = PH.get_node_sets(Dict("Discretization" => raw), dir, mesh)
    @test keys(typed) == keys(reference)
    for key in keys(reference)
        @test typed[key] == reference[key]
    end
end

@testset "external topology file" begin
    dir = mktempdir()
    write(joinpath(dir, "topo.txt"), "1 2 3\n")
    ctx = PS.ParseContext()
    d = PS.convert_value(ID.DiscretizationParams,
                         Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                                          "Input External Topology" => Dict{String,Any}("File" => "topo.txt")),
                         "D", ctx)
    @test PH.external_topology_file(d, dir) == "topo.txt"
    @test_throws PeriLab.PeriLabError PH.external_topology_file(d, mktempdir())
    d2 = PS.convert_value(ID.DiscretizationParams,
                          Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt"),
                          "D", ctx)
    @test PH.external_topology_file(d2, dir) === nothing
end
