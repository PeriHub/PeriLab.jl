import .Exodus
import PeriLab


@timeit to "PeriLab" begin
    PeriLab.main(["strain_xx.yaml"], to)
end

#exodiff("strain_xx.e", "./Reference/ref_strain_xx.e")
#@timeit to "PeriLab" begin
#    PeriLab.main(["strain_xx.yaml"], to)
#end
#
#@timeit to "PeriLab" begin
#    PeriLab.main(["strain_xx.yaml"], to)
#end

