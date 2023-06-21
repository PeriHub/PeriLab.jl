using MPI
using YAML
include("./IO/IO.jl")
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

main()