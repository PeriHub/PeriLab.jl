# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const PH = PeriLab.Parameter_Handling

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

const UT_BM_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const UT_BM_SKIP = ("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml",)

function ut_bm_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_BM_ROOT, root))
        append!(decks, [joinpath(dir, f) for f in files if endswith(f, ".yaml")])
    end
    return sort!(decks)
end

function ut_bm_compare_all_decks()
    for file in ut_bm_decks()
        relpath(file, UT_BM_ROOT) in UT_BM_SKIP && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        deck = raw["PeriLab"]
        input, ctx = ID.read_input(deck, dirname(file))
        deck_needs_license(ctx) && continue
        blocks = input.sections.blocks
        ids = Int64[b["Block ID"] for b in values(deck["Blocks"])]
        @test ID.block_names_and_ids(blocks, ids, true) ==
              PH.get_block_names_and_ids(deck, ids, true)
        @test ID.mesh_scaling(input.sections.discretization) == PH.get_mesh_scaling(deck)
        for id in unique(ids)
            name, block = ID.block_by_id(blocks, id)
            @test block.density == PH.get_density(deck, id)
            @test block.horizon == PH.get_horizon(deck, id)
            @test something(block.fem, false) == PH.get_fem_block(deck, id)
            if block.specific_heat_capacity !== nothing
                @test block.specific_heat_capacity == PH.get_heat_capacity(deck, id)
            end
            for dof in (2, 3)
                if block.angle_x === nothing
                    @test PH.get_angles(deck, id, dof) === nothing
                elseif dof == 2 || (block.angle_y !== nothing && block.angle_z !== nothing)
                    @test ID.block_angles(name, block, dof) == PH.get_angles(deck, id, dof)
                end
            end
        end
    end
end

@testset "block and mesh helpers equal the Dict getters on every shipped deck" begin
    ut_bm_compare_all_decks()
end
