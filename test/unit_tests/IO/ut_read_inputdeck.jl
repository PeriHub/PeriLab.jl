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
    dict = Read_Input_Deck.read_input(filename)
    @test haskey(dict, "data")
    @test dict["data"] == 1
    rm(filename)
    fid = open(filename, "w")
    close(fid)
    @test isnothing(Read_Input_Deck.read_input(filename))
    rm(filename)
end

@testset "ut_read_input_file" begin
    dict = Read_Input_Deck.read_input_file("filename")
    @test isnothing(dict)
    filename = "test.xml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " data: 1")
    close(fid)
    dict = Read_Input_Deck.read_input_file(filename)
    @test isnothing(dict)
    rm(filename)
    filename = "test.yaml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " data: 1")
    close(fid)
    dict = Read_Input_Deck.read_input_file(filename)
    @test isnothing(dict)
    rm(filename)
    filename = "test.yaml"
    fid = open(filename, "w")
    println(fid, "PeriLab:")
    println(fid, " Physics:")
    println(fid, "  d: 3")
    println(fid, "  a: 1")
    println(fid, " Discretization: 1")
    println(fid, " Blocks: 1.7")
    println(fid, " Solver: true")
    close(fid)
    dict = Read_Input_Deck.read_input_file(filename)
    @test dict["Physics"]["d"] == 3
    @test dict["Physics"]["a"] == 1
    @test dict["Discretization"] == 1
    @test dict["Blocks"] == 1.7
    @test dict["Solver"]
    rm(filename)
end
