# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus
using Test
using MPI

export run_perilab

function run_perilab(filename, cores, compare, folder_name="")
    main_path = dirname(@__FILE__)[1:end-4] * "src/main.jl"
    command = `$(Base.julia_cmd()) $main_path -s $(filename).yaml`
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