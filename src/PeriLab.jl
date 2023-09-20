"Simple dummy project module to demonstrate how a project can be organized."
module PeriLab
# using PrecompileTools
# @compile_workload begin
include("./Support/data_manager.jl")
include("./IO/IO.jl")
include("./Core/Solver/Solver_control.jl")
using MPI
using TimerOutputs
using LoggingExtras
const to = TimerOutput()
# using HDF5
using .Data_manager
import .IO
import .Solver
# end
using ArgParse


export main

function print_banner()
    println("""

    8888888b.                  d8b 888               888       |  Documentation: https://docs.julialang.org
    888   Y88b                 Y8P 888               888       |
    888    888                     888               888       |  Type "?" for help, "]?" for Pkg help.
    888   d88P .d88b.  888d888 888 888       8888b.  88888b.   |
    8888888P" d8P  Y8b 888P"   888 888          "88b 888 "88b  |  Version 1.9.1 (2023-06-07)
    888       88888888 888     888 888      .d888888 888  888  |
    888       Y8b.     888     888 888      888  888 888 d88P  |  Official https://julialang.org/ release
    888        "Y8888  888     888 88888888 "Y888888 88888P"   |                                                  

    Copyright: Dr.-Ing. Christian Willberg
    Contact:   christian.willberg@dlr.de

    """)
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
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
        println("Parsed args:")
        for (arg, val) in parsed_args
            println("  $arg  =>  $val")
        end
    end
    main(parsed_args["filename"], parsed_args["dry_run"], parsed_args["verbose"], parsed_args["debug"], parsed_args["silent"])
end
"""
    main(filename, to, dry_run, verbose)

Main function that performs the core functionality of the program.

# Arguments
- `filename`: The name of the file to process.
- `to`: The destination directory.
- `dry_run`: If `true`, performs a dry run without actually moving the file. Default is `false`.
- `verbose`: If `true`, prints additional information during the execution. Default is `false`.
"""
function main(filename, dry_run=false, verbose=false, debug=false, silent=false)

    if debug
        demux_logger = TeeLogger(
            MinLevelLogger(FileLogger(split(filename, ".")[1] * ".log"), Logging.Debug),
            MinLevelLogger(ConsoleLogger(stderr), Logging.Debug),
        )
    else
        demux_logger = TeeLogger(
            MinLevelLogger(FileLogger(split(filename, ".")[1] * ".log"), Logging.Info),
            MinLevelLogger(ConsoleLogger(stderr), Logging.Info),
        )
    end
    global_logger(demux_logger)

    @timeit to "PeriLab" begin
        # init MPI as always ...
        MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)
        if rank == 0 && !silent
            print_banner()
        else
            Logging.disable_logging(Logging.Error)
        end
        #global juliaPath = Base.Filesystem.pwd() * "/"
        global juliaPath = "./"

        ################################
        filename = juliaPath * filename
        # @info filename

        @timeit to "IO.initialize_data" datamanager, params = IO.initialize_data(filename, Data_manager, comm, to)

        @timeit to "Solver.init" blockNodes, bcs, datamanager, solver_options = Solver.init(params, datamanager)
        @timeit to "IO.init_write_results" exos, outputs = IO.init_write_results(params, datamanager)

        # h5write("/tmp/test.h5", "solver_options", solver_options)
        # h5write("/tmp/test.h5", "blockNodes", blockNodes)
        # h5write("/tmp/test.h5", "bcs", 1)
        # h5write("/tmp/test.h5", "datamanager", datamanager)
        # h5write("/tmp/test.h5", "outputs", outputs)
        # h5write("/tmp/test.h5", "exos", exos)

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
        MPI.Finalize()

        if rank == 0
            IO.merge_exodus_files(exos)
        end
        if size > 1
            # IO.delete_files(exos)
        end
    end

    if verbose
        @info demux_logger to
    end
end

# main()

end # module