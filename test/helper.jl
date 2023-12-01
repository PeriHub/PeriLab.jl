# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus
using Test
using MPI
using JSON3
using Logging

export run_perilab
export run_mpi_test
export push_test!

function run_perilab(filename, cores, compare, folder_name="")
    main_path = dirname(@__FILE__)[1:end-4] * "src/main.jl"
    image_path = dirname(@__FILE__)[1:end-4] * "PeriLab_Image"
    command = `$image_path -s $(filename).yaml`
    if cores == 1
        exit_code = run(command).exitcode
        @test exit_code == 0
    else
        mpiexec() do exe  # MPI wrapper

            cmd = `$exe -n $cores $command`
            exit_code = run(cmd).exitcode
            @test exit_code == 0
        end
    end
    if compare
        same = exodiff(filename * ".e", "./Reference/" * filename * ".e"; command_file=folder_name * ".cmd")
        @test same
        if same
            rm("exodiff.log")
            rm(filename * ".e")
            # rm(filename * ".log")
        else
            mv("exodiff.log", filename * "_exodiff.log", force=true)
        end
    end
end

function run_mpi_test(filename, cores, check, folder_name="")
    Logging.disable_logging(Logging.Error)
    command = `$(Base.julia_cmd()) $(folder_name)/$(filename)`
    if cores == 1
        exit_code = run(command).exitcode
        @test exit_code == 0
    else
        mpiexec() do exe  # MPI wrapper

            cmd = `$exe -n $cores $command`
            exit_code = run(cmd).exitcode
            @test exit_code == 0
        end
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

function push_test!(dict::Dict, test::Bool, file::String, line::Int)
    push!(dict["tests"], test)
    push!(dict["line"], "$file:$line")
end