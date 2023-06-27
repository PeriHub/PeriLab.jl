"Simple dummy project module to demonstrate how a project can be organized."
module PeriLab

import MPI
using YAML
include("./IO/IO.jl")

export main
"""
    Main
"""
function main()
    # init MPI as always ...

    MPI.Init()
    global juliaPath = Base.Filesystem.pwd() * "/"
    # from outside #################
    filename::String = "JohnsonCook.yaml"
    comm = MPI.COMM_WORLD

    ################################
    filename = juliaPath * filename
    parameter = read_input_file(filename)
    load_mesh_and_distribute(parameter, comm)
    #

    #include(solver)
    #

    #include("parallel_nn.jl")
    return MPI.Finalize()
end

end # module