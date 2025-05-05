# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#using PeriLab
# using ArgParse

# @testset "ut_parse_commandline" begin
#     @test_throws ArgParseError PeriLab.parse_commandline()
# end

@testset "ut_print_banner" begin
    PeriLab.print_banner(false)
end

@testset "ut_get_examples" begin
    PeriLab.get_examples()
    @test isdir("examples")
    @test isdir("examples/DCB")
    @test isfile("examples/DCB/DCBmodel.yaml")
    rm("examples", recursive = true)
end
