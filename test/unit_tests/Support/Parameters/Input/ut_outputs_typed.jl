# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
function ut_ot(T, raw)
    ctx = PS.ParseContext()
    value = PS.convert_value(T, raw, "X", ctx)
    @test isempty(ctx.errors)
    return value
end

ut_ot_outputs(raw) = ut_ot(Dict{String,ID.OutputParams}, raw)
ut_ot_vars() = Dict{String,Any}("Displacements" => true)

@testset "output frequencies" begin
    outputs = ut_ot_outputs(Dict{String,Any}("freq" => Dict{String,Any}("Output Filename" => "a",
                                                                         "Output Frequency" => 5,
                                                                         "Output Variables" => ut_ot_vars()),
                                             "steps" => Dict{String,Any}("Output Filename" => "b",
                                                                          "Number of Output Steps" => 3,
                                                                          "Output Variables" => ut_ot_vars()),
                                             "both" => Dict{String,Any}("Output Filename" => "c",
                                                                         "Output Frequency" => 1,
                                                                         "Number of Output Steps" => 2,
                                                                         "Output Variables" => ut_ot_vars()),
                                             "per_step" => Dict{String,Any}("Output Filename" => "d",
                                                                             "Output Frequency" => "10 100",
                                                                             "Output Variables" => ut_ot_vars()),
                                             "clamped" => Dict{String,Any}("Output Filename" => "e",
                                                                            "Output Frequency" => 50,
                                                                            "Output Variables" => ut_ot_vars())))
    by_name(nsteps, step) = Dict(zip(keys(outputs),
                                     ID.output_frequencies(outputs, nsteps, step)))
    f = by_name(20, 1)
    @test f["freq"] == 5
    @test f["steps"] == 7          # ceil(20 / 3)
    @test f["both"] == 10          # Number of Output Steps wins: ceil(20 / 2)
    @test f["per_step"] == 10      # first entry for step 1
    @test f["clamped"] == 20       # clamped to nsteps
    @test by_name(200, 2)["per_step"] == 100
end

@testset "filenames and frequencies follow the same output order" begin
    outputs = ut_ot_outputs(Dict{String,Any}("o$i" => Dict{String,Any}("Output Filename" => "file_$i",
                                                                       "Output Frequency" => i,
                                                                       "Output File Type" => isodd(i) ?
                                                                                            "Exodus" :
                                                                                            "CSV",
                                                                       "Output Variables" => ut_ot_vars())
                                             for i in 1:6))
    filenames = ID.output_filenames(outputs, "out")
    frequencies = ID.output_frequencies(outputs, 100, 1)
    for (k, output) in enumerate(values(outputs))
        expected = joinpath("out",
                            output.output_filename *
                            (output.output_file_type == "CSV" ? ".csv" : ".e"))
        @test filenames[k] == expected
        @test frequencies[k] == output.output_frequency
    end
    duplicate = ut_ot_outputs(Dict{String,Any}("a" => Dict{String,Any}("Output Filename" => "same",
                                                                       "Output Frequency" => 1,
                                                                       "Output Variables" => ut_ot_vars()),
                                               "b" => Dict{String,Any}("Output Filename" => "same",
                                                                       "Output Frequency" => 1,
                                                                       "Output Variables" => ut_ot_vars())))
    @test_throws PeriLab.PeriLabError ID.output_filenames(duplicate, "out")
end

@testset "output fieldnames" begin
    variables = Dict("Displacements" => true, "Forces" => true, "Temperature" => false,
                     "Missing" => true, "Reaction" => true)
    field_keys = ["DisplacementsNP1", "Forces"]
    fields = ID.output_fieldnames(variables, field_keys, ["Reaction"], "Exodus")
    @test sort(fields) == sort([["Displacements", "NP1"], ["Forces", "Constant"],
                                ["Reaction", "Constant"]])
    @test ID.output_fieldnames(variables, field_keys, ["Reaction"], "CSV") ==
          [["Reaction", "Constant"]]
end

@testset "active computes" begin
    computes = ut_ot(Dict{String,ID.ComputeClassParams},
                     Dict{String,Any}("b_max" => Dict{String,Any}("Compute Class" => "Block_Data",
                                                                  "Variable" => "Displacements",
                                                                  "Block" => "block_1"),
                                      "a_sum" => Dict{String,Any}("Compute Class" => "Block_Data",
                                                                  "Variable" => "Forces",
                                                                  "Block" => "block_1"),
                                      "bad" => Dict{String,Any}("Compute Class" => "Block_Data",
                                                                "Variable" => "Nope",
                                                                "Block" => "block_1")))
    @test ID.compute_names(computes) == ["a_sum", "b_max", "bad"]
    active = ID.active_computes(computes, ["DisplacementsNP1", "Forces"])
    @test sort(collect(keys(active))) == ["a_sum", "b_max"]
    @test active["b_max"] === computes["b_max"]
end

@testset "ut_check_for_duplicates" begin
    @test !(PeriLab.InputDeck.check_for_duplicates(["a", "b", "c"]))
    @test_logs (:error, "Filename a is used 2 times") @test_throws PeriLab.PeriLabError begin
        PeriLab.InputDeck.check_for_duplicates([
                                                            "a",
                                                            "b",
                                                            "c",
                                                            "a"
                                                        ])
    end
end

