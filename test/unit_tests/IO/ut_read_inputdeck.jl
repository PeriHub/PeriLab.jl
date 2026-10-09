# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
#using PeriLab
@testset "ut_read_input" begin
    filename = "test.yaml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " data: 1")
    close(fid)
    dict = PeriLab.IO.read_input(filename)
    @test haskey(dict["PeriLab"], "data")
    @test dict["PeriLab"]["data"] == 1
    rm(filename)

    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, "data")
    close(fid)
    @test_logs (:error, "Failed to read $filename.") @test_throws PeriLab.PeriLabError begin
        PeriLab.IO.read_input(filename)
    end
    rm(filename)
end

@testset "ut_read_input_deck" begin
    @test_logs (:error,
                "filename can not be found. Make sure the file exist and is readable.") @test_throws PeriLab.PeriLabError begin
        PeriLab.IO.read_input_deck("filename")
    end
    filename = "test.xml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " data: 1")
    close(fid)
    @test_logs (:error,
                "Not a supported filetype $filename") @test_throws PeriLab.PeriLabError begin
        PeriLab.IO.read_input_deck(filename)
    end
    rm(filename)
    filename = "test.yaml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " Models:")
    println(fid, "  d: 3")
    println(fid, "  a: 1")
    println(fid, " Discretization:")
    println(fid, "  Input Mesh File: test")
    println(fid, "  Type: test")
    println(fid, " Blocks:")
    println(fid, "  Block_1:")
    println(fid, "   Block ID: 1")
    println(fid, "   Density: 1.0")
    println(fid, "   Horizon: 1.0")
    println(fid, " Solver:")
    println(fid, "  Initial Time: 0.0")
    println(fid, "  Final Time: 1.0")
    println(fid, "  Verlet:")
    println(fid, "   Safety Factor: 1.0")
    close(fid)
    input = PeriLab.IO.read_input_deck(filename; no_strict = true)  # Models keys d, a are placeholders
    @test input isa PeriLab.InputDeck.PeriLabInput
    @test input.sections.discretization.input_mesh_file == "test"
    @test input.sections.discretization.type == "test"
    @test input.sections.solver.initial_time == 0.0
    @test input.sections.solver.final_time == 1.0
    rm(filename)
end

@testset "read_input_deck errors" begin
    @test_logs (:error,
                "missing.yaml can not be found. Make sure the file exist and is readable.") @test_throws PeriLab.PeriLabError PeriLab.IO.read_input_deck("missing.yaml")
    file = joinpath(mktempdir(), "deck.txt")
    write(file, "x")
    @test_logs (:error, "Not a supported filetype $file") @test_throws PeriLab.PeriLabError PeriLab.IO.read_input_deck(file)
end
