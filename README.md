<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# PeriLab
PeriLab is a Software to solve Peridynamic problems. It is written in Julia to overcome many issues related to the Software Peridigm.
## Modules
```
start julia in batch
in julia promp write ]
now you can define your enviroment via "activate ."
write instantiate

to add the modules you can use 

add MPI
add NearestNeighbors
add Logging
add YAML
add CSV
add DataFrames
add Exodus
add LinearAlgebra


```

```
If you want to define it in an initialise file use

import Pkg
Pkg.add("MPI")
Pkg.add("NearestNeighbors")
Pkg.add("Logging")
Pkg.add("YAML")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Exodus")
Pkg.add("LinearAlgebra")
```
## Run PeriLab
 "/home/$user/.julia/bin/mpiexecjl  -n $ncores /home/$user/julia/julia-1.9.1/bin/julia ./src/main.jl""sr

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