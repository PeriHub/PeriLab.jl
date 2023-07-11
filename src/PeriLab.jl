"Simple dummy project module to demonstrate how a project can be organized."
module PeriLab
include("./Support/data_manager.jl")
include("./IO/IO.jl")
import MPI
using YAML
import .Data_manager


export main
"""
    Main
"""
function main()
    # init MPI as always ...

    MPI.Init()
    #global juliaPath = Base.Filesystem.pwd() * "/"
    global juliaPath = "./"
    # from outside #################
    filename::String = "Input.yaml"
    comm = MPI.COMM_WORLD

    ################################
    filename = juliaPath * filename

    datamanager, parameter = init_data(filename, Data_manager, comm)


    #    save_values(results)
    return MPI.Finalize()
end

end # module