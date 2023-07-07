using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
    Pkg.instantiate()
end
import PeriLab

PeriLab.main()