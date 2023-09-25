using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")

    Pkg.instantiate()
end
# using Revise
import PeriLab
global PeriLabPath = "/home/will_cr/julia_projects/perilab/"
PeriLab.main()
