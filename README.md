<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# PeriLab
PeriLab is a Software to solve Peridynamic problems. It is written in Julia to overcome many issues related to the Software Peridigm.

## Documentation

https://dlr-perihub.gitlab.io/perilab/

## Installation

The PeriLab package is available through the Julia package system and can be installed using the following commands:

```julia
using Pkg
Pkg.add("PeriLab")
```

Throughout the rest of this tutorial, we will assume that you have installed the
PeriLab package and have already typed `using PeriLab` to bring all of the
relevant variables into your current namespace.

## Using `PeriLab` 

The simplest way to run the `PeriLab` simulation core is to use the provided examples. 

```julia PeriLab
using PeriLab

PeriLab.get_examples()
PeriLab.main("examples/DCB/DCBmodel.yaml")
```
The main functionalities for the `yaml` input deck is given in
```
"examples/functionalities.yaml"
```

## Using `PeriLab` with multiple processors (MPI)

In order to run `PeriLab` for large scale problems [MPI](https://juliaparallel.org/MPI.jl/stable/usage/) needs to be installed:

```sh
$ julia
julia> using MPI
julia> MPI.install_mpiexecjl()
```

Run PeriLab with two processors:
```sh
$ mpiexecjl -n 2 julia -e "using PeriLab; PeriLab.main()" examples/DCB/DCBmodel.yaml -v
```

## Contributing

Contributions in all forms (bug reports, documentation, features, suggestions, ...) are very
welcome. 

## Questions
If you have questions about PeriLab.jl you're welcome to reach out to us on the authors e mail adresses.

## Authors and acknowledgment
[Dr.-Ing. Christian Willberg](mailto::christian.willberg@dlr.de)

[M.Sc. Jan-Timo Hesse](mailto::jan-timo.hesse@dlr.de)

## Project status
In development.
