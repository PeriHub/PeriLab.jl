"Simple dummy project module to demonstrate how a project can be organized."
module PeriLab

import MPI
using YAML
import .Parameter_Handling
include("./IO/IO.jl")

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

    datamanager, parameter = init_data(filename, Parameter_Handling, comm)
    #results = solver(data, parameter, comm)

    #    save_values(results)
    return MPI.Finalize()
end

end # module