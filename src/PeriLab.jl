"Simple dummy project module to demonstrate how a project can be organized."
module PeriLab
include("./Support/data_manager.jl")
include("./IO/IO.jl")
include("./Core/Solver/Solver_control.jl")
using MPI
using LoggingExtras
using .Data_manager
import .IO
import .Solver

export main
"""
    Main
"""
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
function main(ARGS)
    demux_logger = TeeLogger(
        MinLevelLogger(FileLogger(split(ARGS[1], ".")[1] * ".log"), Logging.Info),
        MinLevelLogger(ConsoleLogger(stderr), Logging.Info),
    )
    global_logger(demux_logger)
    # init MPI as always ...

    MPI.Init()
    comm = MPI.COMM_WORLD
    if MPI.Comm_rank(comm) == 0
        print_banner()
    end
    #global juliaPath = Base.Filesystem.pwd() * "/"
    #global juliaPath = "./"
    # from outside #################

    @info ARGS
    filename::String = ARGS[1]
    # filename::String = "Input.yaml"

    ################################
    #filename =  filename
    # filename = "Input.yaml"
    datamanager, params = IO.initialize_data(filename, Data_manager, comm)

    blockNodes, bcs, datamanager, solver_options = Solver.init(params, datamanager)
    exos, outputs = IO.init_write_results(params, datamanager)
    exos = Solver.solver(solver_options, blockNodes, bcs, datamanager, outputs, exos, IO.write_results)

    IO.close_files(exos)

    return MPI.Finalize()
end

end # module