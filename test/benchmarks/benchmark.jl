using TimerOutputs
using JSON3
const to = TimerOutput()

include("../helper.jl")

cd("test/benchmarks") do
    @timeit to "perilab 1" run_perilab("strain_xx", 1, false)
    @timeit to "perilab 2" run_perilab("strain_xx", 2, false)
    @timeit to "perilab 3" run_perilab("strain_xx", 3, false)
    @timeit to "perilab 4" run_perilab("strain_xx", 4, false)
end

open("benchmark_results.json", "w") do f
    JSON3.pretty(f, JSON3.write(TimerOutputs.todict(to)))
end