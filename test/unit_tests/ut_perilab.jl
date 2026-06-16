# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
#using PeriLab
using Logging
# using ArgParse: ArgParseError

@testset "ut_print_banner" begin
    PeriLab.print_banner(false)
end

@testset "ut_parse_commandline" begin
    @test PeriLab.parse_commandline([]) ==
          Dict("dry_run" => false, "reload" => false, "examples" => false,
               "verbose" => false, "output_dir" => "", "filenames" => [], "silent" => false,
               "debug" => false)
    @test PeriLab.parse_commandline(["Test.yaml", "-r", "-v", "-d", "-s", "-o", "folder"]) ==
          Dict("dry_run" => false, "reload" => true, "examples" => false, "verbose" => true,
               "output_dir" => "folder", "filenames" => ["Test.yaml"], "silent" => true,
               "debug" => true)
end

@testset "ut_main" begin
    @test_logs (:error, "Please provide at least one filename, f.e. 'PeriLab example.yaml'") PeriLab.main([])
    @test_logs min_level=Logging.Warn PeriLab.main([
                                                       "fullscale_tests/test_bond_based_elastic/test_bond_based_elastic_1D.yaml",
                                                       "-s"
                                                   ])
end

# @testset "ut_get_examples" begin
#     PeriLab.get_examples()
#     @test isdir("examples")
#     @test isdir("examples/DCB")
#     @test isfile("examples/DCB/DCBmodel.yaml")
#     rm("examples", recursive = true)
# end
