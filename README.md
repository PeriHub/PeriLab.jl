<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# PeriLab
PeriLab is a Software to solve Peridynamic problems. It is written in Julia to overcome many issues related to the Software Peridigm.

## Documentation

https://fa_sw.pages.gitlab.dlr.de/peridynamik/perilab/

## Installation
```
download repository
go in existing_repo
start julia in batch
in julia promp write ]
now you can define your enviroment via "activate ."
write instantiate
```

## Run PeriLab

julia --project=../../ ../../src/main.jl dcb.yaml -v

or for MPI

 "/home/$user/.julia/bin/mpiexecjl  -n $ncores /home/$user/julia/julia-1.9.1/bin/julia ./src/main.jl""sr

## Contributing

Contributions in all forms (bug reports, documentation, features, suggestions, ...) are very
welcome. 

## Questions
If you have questions about Ferrite.jl you're welcome to reach out to us on the authors e mail adresses.
## Authors and acknowledgment
'''
Dr.-Ing. Christian Willberg; christian.willberg@dlr.de
M.Sc. Jan-Timo Hesse; jan-timo.hesse@dlr.de
'''
## Project status
In development.
