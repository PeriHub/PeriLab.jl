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

