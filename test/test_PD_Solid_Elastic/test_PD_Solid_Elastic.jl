using Exodus
using Test
using MPI
# import PeriLab

function check_exodiff(filename, cores)
    mpiexec() do exe  # MPI wrapper

        run(`$exe -n $cores $(Base.julia_cmd()) ../../src/main.jl -s $(filename).yaml`)
        # cmd = """$exe -n $n $($(Base.julia_cmd())) [...]/main.jl $filename.yaml"""
        # run(cmd)
        # alternatively:
        # p = run(ignorestatus(`...`))
        # @test success(p)
    end
    # PeriLab.main(filename * ".yaml", true)
    exodiff(filename * ".e", "./Reference/" * filename * ".e")
    # exodiff(filename * ".e.2.0", filename * ".e.2.0")
    @test occursin("Files are the same", read("exodiff.log", String))
    rm("exodiff.log")
    rm(filename * ".e")
end

cd("test/test_PD_Solid_Elastic") do
    check_exodiff("strain_xx", 1)
    check_exodiff("strain_xy", 1)
    check_exodiff("strain_yy", 1)
end