"Simple dummy project module to demonstrate how a project can be organized."
module PeriLab
include("./Support/data_manager.jl")
include("./IO/IO.jl")
import MPI
import .Data_manager
include("./Core/Solver/Solver_control.jl")
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
    """)
end
function main()
    # init MPI as always ...

    MPI.Init()
    comm = MPI.COMM_WORLD
    if MPI.Comm_rank(comm) == 0
        print_banner()
    end
    #global juliaPath = Base.Filesystem.pwd() * "/"
    global juliaPath = "./"
    # from outside #################
    filename::String = "Input.yaml"


    ################################
    filename = juliaPath * filename

    datamanager, parameter = IO.initialize_data(filename, Data_manager, comm)


    #write_results(datamanager)


    #    save_values(results)
    return MPI.Finalize()
end

end # module