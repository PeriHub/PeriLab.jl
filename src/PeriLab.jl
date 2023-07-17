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
function main()
    include("banner.jl")
    # init MPI as always ...

    MPI.Init()
    #global juliaPath = Base.Filesystem.pwd() * "/"
    global juliaPath = "./"
    # from outside #################
    filename::String = "Input.yaml"
    comm = MPI.COMM_WORLD

    ################################
    filename = juliaPath * filename

    datamanager, parameter = IO.init_data(filename, Data_manager, comm)

    #write_results(datamanager)


    #    save_values(results)
    return MPI.Finalize()
end

end # module