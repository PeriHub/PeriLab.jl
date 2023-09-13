using Exodus
using Test
using TimerOutputs
import PeriLab

const to = TimerOutput()

cd("test_PD_Solid_Elastic/") do
    PeriLab.main("strain_xx.yaml", to, false, false, true)
    exodiff("strain_xx.e", "./Reference/strain_xx.e")
    @test occursin("Files are the same", read("exodiff.log", String))
    rm("exodiff.log")
    rm("strain_xx.e")
end
#@timeit to "PeriLab" begin
#    PeriLab.main(["strain_xx.yaml"], to)
#end
#
#@timeit to "PeriLab" begin
#    PeriLab.main(["strain_xx.yaml"], to)
#end

