# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus
using Test
using MPI
using JSON3
using Logging
using CSV
using DataFrames
using PeriLab

export run_perilab
export run_and_compare
export run_mpi_test
export push_test!

function run_perilab(filename, cores, compare, folder_name = ""; silent = true,
                     reload = false, compare_csv = false)
    if cores == 1
        PeriLab.main(filename * ".yaml"; silent = silent, reload = reload)
    else
        fn = """using PeriLab; PeriLab.main(""" *
             '"' *
             filename *
             """.yaml"; silent = true)"""
        command = `$(Base.julia_cmd()) -e "$(fn)"`
        mpiexec() do exe  # MPI wrapper
            cmd = `$exe -n $cores $command`
            exit_code = run(cmd).exitcode
            @test exit_code == 0
        end
    end
    if compare
        same = false
        if cores == 1
            same = exodiff(filename * ".e",
                           "./Reference/" * filename * ".e";
                           command_file = folder_name * ".cmd",)
        else
            same = exodiff(filename * ".e",
                           "./Reference/" * filename * ".e",
                           ["-p", "-f", folder_name * ".cmd"])
        end
        @test same
        if same
            rm("exodiff.log")
            rm(filename * ".e")
        else
            mv("exodiff.log", filename * "_exodiff.log", force = true)
        end
    else
        rm(filename * ".e")
    end
    if compare_csv
        results_csv = CSV.read(open(filename * ".csv"), DataFrame)
        reference_csv = CSV.read(open("./Reference/" * filename * ".csv"), DataFrame)
        for i in 1:size(reference_csv)[2]
            @test results_csv[:, i] ≈ reference_csv[:, i] atol=1e-15
        end
        rm(filename * ".csv")
    end
end

function run_and_compare(filename1, filename2)
    PeriLab.main(filename1 * ".yaml"; silent = true)
    PeriLab.main(filename2 * ".yaml"; silent = true)
    exo1 = ExodusDatabase(filename1 * ".e", "r")
    exo2 = ExodusDatabase(filename2 * ".e", "r")

    n_steps1 = read_number_of_time_steps(exo1)
    displ_x_1 = read_values(exo1, NodalVariable, n_steps1, 1, "Displacementsx")
    displ_y_1 = read_values(exo1, NodalVariable, n_steps1, 1, "Displacementsy")
    force_x_1 = read_values(exo1, NodalVariable, n_steps1, 1, "Forcesx")
    force_y_1 = read_values(exo1, NodalVariable, n_steps1, 1, "Forcesy")

    n_steps2 = read_number_of_time_steps(exo2)
    displ_x_2 = read_values(exo2, NodalVariable, n_steps2, 1, "Displacementsx")
    displ_y_2 = read_values(exo2, NodalVariable, n_steps2, 1, "Displacementsy")
    force_x_2 = read_values(exo2, NodalVariable, n_steps2, 1, "Forcesx")
    force_y_2 = read_values(exo2, NodalVariable, n_steps2, 1, "Forcesy")

    @test n_steps1 == n_steps2
    @test displ_x_1 == displ_x_2
    @test displ_y_1 == displ_y_2
    @test force_x_1 == force_x_2
    @test force_y_1 == force_y_2

    if n_steps1 == n_steps2 &&
       displ_x_1 == displ_x_2 &&
       displ_y_1 == displ_y_2 &&
       force_x_1 == force_x_2 &&
       force_y_1 == force_y_2
        rm(filename1 * ".e")
        rm(filename2 * ".e")
    end
end

function run_mpi_test(filename, cores, check, folder_name = "")
    command = `$(Base.julia_cmd()) $(folder_name)/$(filename)`

    mpiexec() do exe  # MPI wrapper
        cmd = `$exe -n $cores $command`
        exit_code = run(cmd).exitcode
        @test exit_code == 0
    end
    if check
        check_test_json(cores)
    end
end

function check_test_json(cores)
    for i in 1:cores
        @testset "core $i" begin
            failed = false
            data = JSON3.read(read("test_results_$(i-1).json", String))
            for testset in keys(data)
                @testset "$testset" begin
                    for test in data[testset]["tests"]
                        @test test
                        if !test
                            failed = true
                        end
                    end
                end
            end
            if !failed
                rm("test_results_$(i-1).json")
            end
        end
    end
end
