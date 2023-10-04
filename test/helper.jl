# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus
using Test
using MPI

export run_perilab

function run_perilab(filename, cores, compare)
    if cores == 1
        cmd = `$(Base.julia_cmd()) ../../src/main.jl -s $(filename).yaml`
        exit_code = run(cmd).exitcode
        @test exit_code == 0
    else
        mpiexec() do exe  # MPI wrapper

            cmd = `$exe -n $cores $(Base.julia_cmd()) ../../src/main.jl -s $(filename).yaml`
            exit_code = run(cmd).exitcode
            @test exit_code == 0
        end
    end
    if compare
        exodiff(filename * ".e", "./Reference/" * filename * ".e")
        @test occursin("Files are the same", read("exodiff.log", String))
        rm("exodiff.log")
    end
    rm(filename * ".e")
    rm(filename * ".log")
end