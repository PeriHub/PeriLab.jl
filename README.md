# PeriLab
PeriLab is a Software to solve Peridynamic problems. It is written in Julia to overcome many issues related to the Software Peridigm.
## Modules
```
start julia in batch
in julia promp write ]
now you can define your enviroment via "activate ."
write instantiate
```

```
import Pkg
Pkg.add("MPI")
Pkg.add("NearestNeighbors")
Pkg.add("Logging")
Pkg.add("YAML")
Pkg.add("CSV")
Pkg.add("DataFrames")
```
## Add your files

```
cd existing_repo
git remote add origin https://gitlab.dlr.de/fa_sw/peridynamik/perilab.git
git branch -M main
git push -uf origin main
```
## Authors and acknowledgment
'''
Dr.-Ing. Christian Willberg; christian.willberg@dlr.de
M.Sc. Jan-Timo Hesse; jan-timo.hesse@dlr.de
'''
## Project status
In development.
