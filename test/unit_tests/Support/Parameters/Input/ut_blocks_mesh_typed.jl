# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
function ut_bm(T, raw)
    ctx = PS.ParseContext()
    value = PS.convert_value(T, raw, "X", ctx)
    return value, ctx
end

ut_bm_blocks(raw) = first(ut_bm(Dict{String,ID.BlockParams}, raw))

const UT_BM_BLOCKS = ut_bm_blocks(Dict{String,Any}("left" => Dict{String,Any}("Block ID" => 1,
                                                                              "Density" => 2.0,
                                                                              "Horizon" => 0.5,
                                                                              "Angle X" => 30.0),
                                                   "right" => Dict{String,Any}("Block ID" => 3,
                                                                               "Density" => 4.0,
                                                                               "Horizon" => 0.7,
                                                                               "Angle X" => 10.0,
                                                                               "Angle Y" => 20.0,
                                                                               "Angle Z" => 0.0)))

@testset "block by id" begin
    name, block = ID.block_by_id(UT_BM_BLOCKS, 3)
    @test name == "right" && block.density === 4.0
    @test_throws PeriLab.PeriLabError ID.block_by_id(UT_BM_BLOCKS, 2)
end

@testset "block angles" begin
    @test ID.block_angles(ID.block_by_id(UT_BM_BLOCKS, 1)..., 2) === 30.0
    @test_throws PeriLab.PeriLabError ID.block_angles(ID.block_by_id(UT_BM_BLOCKS, 1)..., 3)
    @test ID.block_angles(ID.block_by_id(UT_BM_BLOCKS, 3)..., 3) == [10.0, 20.0, 0.0]
    no_angles = ut_bm_blocks(Dict{String,Any}("b" => Dict{String,Any}("Block ID" => 1,
                                                                      "Density" => 1.0,
                                                                      "Horizon" => 1.0)))
    @test ID.block_angles(ID.block_by_id(no_angles, 1)..., 3) === nothing
end

@testset "block names and ids" begin
    names, ids = ID.block_names_and_ids(UT_BM_BLOCKS, [1, 3], true)
    @test names == ["left", "right"] && ids == [1, 3]
    names, ids = ID.block_names_and_ids(UT_BM_BLOCKS, [3], true)   # block 1 not in mesh
    @test names == ["right"] && ids == [3]
end

@testset "mesh scaling" begin
    d, ctx = ut_bm(ID.DiscretizationParams,
                   Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                                    "Horizon Mesh Scaling Y" => 2))
    @test isempty(ctx.errors)
    @test ID.mesh_scaling(d) == [1.0, 2.0, 1.0]
end

@testset "gcode block ids" begin
    gcode(blocks) = Dict{String,Any}("Overwrite Mesh" => true, "Sampling" => 1.0,
                                     "Width" => 1.0, "Height" => 1.0, "Blocks" => blocks)
    g, ctx = ut_bm(ID.GcodeParams, gcode(Dict{Any,Any}(2 => "z > 1", 1 => "z <= 1")))
    @test isempty(ctx.errors)
    @test ID.gcode_block_ids(g) == Dict(2 => "z > 1", 1 => "z <= 1")
    g, ctx = ut_bm(ID.GcodeParams, gcode(Dict{Any,Any}("top" => "z > 1")))
    @test ctx.errors[1].path == "X.Blocks.top"
    @test ctx.errors[1].message == "expected a block id (integer) as key"
    g, _ = ut_bm(ID.GcodeParams, Dict{String,Any}("Overwrite Mesh" => true, "Sampling" => 1.0,
                                                  "Width" => 1.0, "Height" => 1.0))
    @test ID.gcode_block_ids(g) === nothing
end

@testset "bond filter required fields" begin
    _, ctx = ut_bm(ID.BondFilterParams,
                   Dict{String,Any}("Type" => "Disk", "Normal X" => 0.0, "Normal Y" => 0.0,
                                    "Normal Z" => 1.0, "Center X" => 0.0, "Center Y" => 0.0,
                                    "Center Z" => 0.0))
    @test ctx.errors[1].path == "X"
    @test ctx.errors[1].message == "\"Disk\" bond filter requires: Radius"
    _, ctx = ut_bm(ID.BondFilterParams,
                   Dict{String,Any}("Type" => "Rectangular_Plane", "Normal X" => 0.0,
                                    "Normal Y" => 1.0))
    @test ctx.errors[1].message ==
          "\"Rectangular_Plane\" bond filter requires: Lower Left Corner X, Lower Left Corner Y, Bottom Unit Vector X, Bottom Unit Vector Y, Bottom Length, Side Length"
    _, ctx = ut_bm(ID.BondFilterParams,
                   Dict{String,Any}("Type" => "My_Filter", "Normal X" => 0.0, "Normal Y" => 1.0))
    @test isempty(ctx.errors)
end

