using Test
using MPI
using TimerOutputs
using JSON3
const to = TimerOutput()

function run_perilab(filename, cores)
    if cores == 1
        run(`$(Base.julia_cmd()) ../../src/main.jl -v $(filename).yaml`)
    else
        mpiexec() do exe  # MPI wrapper

            run(`$exe -n $cores $(Base.julia_cmd()) ../../src/main.jl -v $(filename).yaml`)
        end
    end
    rm(filename * ".e")
    rm(filename * ".log")
end

cd("test/benchmarks") do
    @timeit to "perilab 1" run_perilab("strain_xx", 1)
    @timeit to "perilab 2" run_perilab("strain_xx", 2)
    @timeit to "perilab 3" run_perilab("strain_xx", 3)
    @timeit to "perilab 4" run_perilab("strain_xx", 4)
end

open("benchmark_results.json", "w") do f
    JSON3.pretty(f, JSON3.write(TimerOutputs.todict(to)))
end