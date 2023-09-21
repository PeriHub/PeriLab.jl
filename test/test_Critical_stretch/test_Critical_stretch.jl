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
    @test occursin("Files are the same", read("exodiff.log", String))
    rm("exodiff.log")
    rm(filename * ".e")
    rm(filename * ".log")
end

cd("test_PD_Solid_Elastic") do
    check_exodiff("critical_stretch_tension", 1)
    check_exodiff("critical_stretch_pressure", 1)
end