include("../helper.jl")

cd(basename(@__FILE__)[1:end-3]) do
    run_perilab("critical_stretch_tension", 1, true)
    run_perilab("critical_stretch_pressure", 1, true)
end