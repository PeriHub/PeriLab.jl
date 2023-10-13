# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"Simple dummy project module to demonstrate how a project can be organized."
module PeriLab
include("./Support/data_manager.jl")
include("./IO/logging.jl")
include("./IO/IO.jl")
include("./Core/Solver/Solver_control.jl")
using MPI
using Pkg
using TimerOutputs
using Logging
const to = TimerOutput()
using .Data_manager
import .Logging_module
import .IO
import .Solver
# end
using ArgParse


export main

function print_banner()
    println("""\e[]
    \e[1;36mPeriLab. \e[0m                  \e[1;32md8b \e[1;36m888               888\e[0m       |  Version: """ * string(Pkg.project().version) * """ 
    \e[1;36m888   Y88b\e[0m                 \e[1;32mY8P \e[1;36m888               888\e[0m       |
    \e[1;36m888    888\e[0m                     \e[1;36m888               888\e[0m       |  Copyright: Dr.-Ing. Christian Willberg, M.Sc. Jan-Timo Hesse 
    \e[1;36m888   d88P\e[0m \e[1;36m.d88b.\e[0m  \e[1;36m888d888 888 888       \e[1;36m8888b.\e[0m  \e[1;36m88888b.\e[0m   |  Contact: christian.willberg@dlr.de, jan-timo.hesse@dlr.de
    \e[1;36m8888888P"\e[0m \e[1;36md8P  Y8b\e[0m \e[1;36m888P"   888 888          \e[1;36m"88b\e[0m \e[1;36m888 "88b\e[0m  |  
    \e[1;36m888\e[0m       \e[1;36m88888888\e[0m \e[1;36m888\e[0m     \e[1;36m888\e[0m \e[1;36m888\e[0m      \e[1;36m.d888888\e[0m \e[1;36m888  888\e[0m  |  License: BSD-3-Clause
    \e[1;36m888\e[0m       \e[1;36mY8b.\e[0m     \e[1;36m888\e[0m     \e[1;36m888\e[0m \e[1;36m888\e[0m      \e[1;36m888  888\e[0m \e[1;36m888 d88P\e[0m  |  
    \e[1;36m888\e[0m        \e[1;36m"Y8888\e[0m  \e[1;36m888\e[0m     \e[1;36m888\e[0m \e[1;36m88888888\e[0m \e[1;36m"Y888888\e[0m \e[1;36m88888P"\e[0m   |  Gitlab: https://gitlab.com/dlr-perihub/perilab                                                
    """)
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--dry_run"
        help = "dry_run"
        action = :store_true
        "--verbose", "-v"
        help = "verbose"
        action = :store_true
        "--debug", "-d"
        help = "debug"
        action = :store_true
        "--silent", "-s"
        help = "silent"
        action = :store_true
        "filename"
        help = "filename"
        required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    if parsed_args["verbose"]
        @info "Parsed args:"
        for (arg, val) in parsed_args
            @info "  $arg  =>  $val"
        end
    end
    main(parsed_args["filename"], parsed_args["dry_run"], parsed_args["verbose"], parsed_args["debug"], parsed_args["silent"])
end

"""
    main(filename, dry_run, verbose, debug, silent)

Main function that performs the core functionality of the program.
# Arguments
- `filename`: The name of the file to process.
- `to`: The destination directory.
- `dry_run`: If `true`, performs a dry run without actually moving the file. Default is `false`.
- `verbose`: If `true`, prints additional information during the execution. Default is `false`.
"""
function main(filename, dry_run=false, verbose=false, debug=false, silent=false)

    @timeit to "PeriLab" begin

        if !silent
            @timeit to "Logging.init_logging" Logging_module.init_logging(filename, debug)
        end

        # init MPI as always ...
        MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)
        if rank == 0 && !silent
            print_banner()
            @info "\n PeriLab version: " * string(Pkg.project().version) * "\n Copyright: Dr.-Ing. Christian Willberg, M.Sc. Jan-Timo Hesse\n Contact: christian.willberg@dlr.de, jan-timo.hesse@dlr.de\n Gitlab: https://gitlab.com/dlr-perihub/perilab\n doi: \n License: BSD-3-Clause\n ---------------------------------------------------------------"
        else
            Logging.disable_logging(Logging.Error)
        end
        #global juliaPath = Base.Filesystem.pwd() * "/"
        global juliaPath = "./"

        ################################
        filename = juliaPath * filename
        # @info filename

        @timeit to "IO.initialize_data" datamanager, params = IO.initialize_data(filename, Data_manager, comm, to)
        @info "Solver init"
        @timeit to "Solver.init" blockNodes, bcs, datamanager, solver_options = Solver.init(params, datamanager)
        @info "Init write results"
        @timeit to "IO.init_write_results" exos, outputs = IO.init_write_results(params, datamanager, solver_options["nsteps"])

        if dry_run
            nsteps = solver_options["nsteps"]
            solver_options["nsteps"] = 10
            elapsed_time = @elapsed begin
                @timeit to "Solver.solver" exos = Solver.solver(solver_options, blockNodes, bcs, datamanager, outputs, exos, IO.write_results, to, silent)
            end

            @info "Estimated runtime: " * string((elapsed_time / 10) * nsteps) * " [s]"
            file_size = IO.get_file_size(exos)
            @info "Estimated filesize: " * string((file_size / 10) * nsteps) * " [b]"

        else
            @timeit to "Solver.solver" exos = Solver.solver(solver_options, blockNodes, bcs, datamanager, outputs, exos, IO.write_results, to, silent)
        end

        IO.close_files(exos)

        if dry_run
            IO.delete_files(exos)
        end

        if size > 1 && rank == 0
            IO.merge_exodus_files(exos)
        end
        MPI.Barrier(comm)
        if size > 1
            IO.delete_files(exos)
        end
        MPI.Finalize()
    end

    if verbose
        @info to
    end
end

end # module