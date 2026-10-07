# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_bc(raw)
    ctx = PS.ParseContext()
    value = PS.convert_value(ID.BoundaryConditionParams, raw, "BC", ctx)
    return value, ctx
end

ut_bc_raw(extra...) = Dict{String,Any}("Variable" => "Displacements", "Node Set" => "Set-1",
                                       "Value" => 0.0, extra...)

@testset "node set names" begin
    bc, _ = ut_bc(ut_bc_raw("Node Set" => "top + left+right"))
    @test ID.bc_node_set_names(bc) == ["top", "left", "right"]
    bc, _ = ut_bc(ut_bc_raw())
    @test ID.bc_node_set_names(bc) == ["Set-1"]
end

@testset "step ids" begin
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw()))) === nothing
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw("Step ID" => 2)))) == [2]
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw("Step ID" => "1,3")))) == [1, 3]
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw("Step ID" => "1, 3")))) == [1, 3]
    _, ctx = ut_bc(ut_bc_raw("Step ID" => "1,a"))
    @test ctx.errors[1].path == "BC.\"Step ID\""
    @test ctx.errors[1].message ==
          "expected an integer or a comma-separated list of integers, got \"1,a\""
end

const UT_BC_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const UT_BC_SKIP = ("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml",)

function ut_bc_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_BC_ROOT, root))
        append!(decks, [joinpath(dir, f) for f in files if endswith(f, ".yaml")])
    end
    return sort!(decks)
end

function ut_bc_compare_all_decks()
    for file in ut_bc_decks()
        relpath(file, UT_BC_ROOT) in UT_BC_SKIP && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        deck = raw["PeriLab"]
        input, ctx = ID.read_input(deck, dirname(file))
        deck_needs_license(ctx) && continue
        reference = PeriLab.Parameter_Handling.get_bc_definitions(deck)
        @test sort(collect(keys(input.sections.boundary_conditions))) ==
              sort(string.(collect(keys(reference))))
        for (name, bc) in input.sections.boundary_conditions
            old = reference[name]
            @test ID.bc_node_set_names(bc) == String.(strip.(split(old["Node Set"], "+")))
            if haskey(old, "Step ID")
                @test ID.bc_step_ids(bc) ==
                      parse.(Int64, strip.(split(string(old["Step ID"]), ",")))
            else
                @test ID.bc_step_ids(bc) === nothing
            end
        end
    end
end

@testset "BC helpers equal the Dict definitions on every shipped deck" begin
    ut_bc_compare_all_decks()
end
