# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/IO/read_inputdeck.jl")

using Test

@testset "ut_read_input" begin
    filename = "test.yaml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " data: 1")
    close(fid)
    dict = read_input(filename)
    @test haskey(dict["PeriLab"], "data")
    @test dict["PeriLab"]["data"] == 1
    rm(filename)
    fid = open(filename, "w")
    close(fid)
    @test isnothing(read_input(filename))
    rm(filename)

    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, "data")
    close(fid)
    @test isnothing(read_input(filename))
    rm(filename)
end

@testset "ut_read_input_file" begin
    dict = read_input_file("filename")
    @test isnothing(dict)
    filename = "test.xml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " data: 1")
    close(fid)
    dict = read_input_file(filename)
    @test isnothing(dict)
    rm(filename)
    filename = "test.yaml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " data: 1")
    close(fid)
    dict = read_input_file(filename)
    @test isnothing(dict)
    rm(filename)
    filename = "test.yaml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " Physics:")
    println(fid, "  d: 3")
    println(fid, "  a: 1")
    println(fid, " Discretization:")
    println(fid, "  Input Mesh File: test")
    println(fid, "  Type: test")
    println(fid, " Blocks:")
    println(fid, "  Block_1:")
    println(fid, " Solver:")
    println(fid, "  Initial Time: 0.0")
    println(fid, "  Final Time: 1.0")
    close(fid)
    dict = read_input_file(filename)
    @test dict["Physics"]["d"] == 3
    @test dict["Physics"]["a"] == 1
    @test dict["Discretization"]["Input Mesh File"] == "test"
    @test dict["Discretization"]["Type"] == "test"
    @test dict["Solver"]["Initial Time"] == 0.0
    @test dict["Solver"]["Final Time"] == 1.0
    rm(filename)
end
