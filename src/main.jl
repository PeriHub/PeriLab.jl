using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")

    Pkg.instantiate()
end
using Revise
import PeriLab
using TimerOutputs
const to = TimerOutput()

@timeit to "PeriLab" begin
    PeriLab.main(ARGS, to)
end

show(to)