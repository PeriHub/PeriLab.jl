using Exodus
using Test
import PeriLab

function check_exodiff(filename)
    PeriLab.main(filename * ".yaml", true)
    exodiff(filename * ".e", "./Reference/" * filename * ".e")
    @test occursin("Files are the same", read("exodiff.log", String))
    rm("exodiff.log")
    rm(filename * ".e")
end

cd("test_PD_Solid_Elastic/") do
    check_exodiff("strain_xx")
    # check_exodiff("strain_xy")
    # check_exodiff("strain_yy")
end