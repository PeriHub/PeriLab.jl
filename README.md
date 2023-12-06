<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# `PeriLab` - Peridynamic Laboratory
Welcome to `PeriLab`, a powerful software solution designed for tackling Peridynamic problems. 

## Documentation

Explore the comprehensive [documentation](https://dlr-perihub.gitlab.io/PeriLab.jl/) for `PeriLab`

## Installation

The `PeriLab`  package is available through the Julia package system and can be installed using the following commands:

```julia
using Pkg
Pkg.add("PeriLab")
```

Throughout the rest of this tutorial, we will assume that you have installed the
PeriLab package and have already typed `using PeriLab` to bring all of the
relevant variables into your current namespace.

## Getting Started with `PeriLab` 

Jumpstart your exploration of the PeriLab simulation core with provided examples. Run the following commands in Julia:

```julia PeriLab
using PeriLab

PeriLab.get_examples()
PeriLab.main("examples/DCB/DCBmodel.yaml")
```
>Note: More details about the main functionalities in the yaml input deck [here](https://gitlab.com/dlr-perihub/PeriLab.jl/-/blob/main/src/Support/Parameters/parameter_handling.jl?ref_type=heads).

## Parallel Processing with `PeriLab` (MPI)

To handle large-scale problems efficiently, install [MPI](https://juliaparallel.org/MPI.jl/stable/usage/). Run PeriLab with two processors using the following commands:

```sh
$ julia
julia> using MPI
julia> MPI.install_mpiexecjl()
```

Run PeriLab with two processors:
```sh
$ mpiexecjl -n 2 julia -e "using PeriLab; PeriLab.main()" examples/DCB/DCBmodel.yaml -v
```

## `PeriLab` on `JuliaHub`

Experience the convenience of using PeriLab as a ready-to-use application on JuliaHub. Simply create an [account](https://juliahub.com), navigate to the [applications page](https://juliahub.com/ui/Applications), and add the repository URL: https://gitlab.com/dlr-perihub/PeriLab.jl.

Configure advanced options, such as _filename_, _dryrun_, _verbosity_, _debug_, and _silence_. Click __Start__ and monitor the job progress. Results will be available in a zipped folder.

Hit the __Start__ button and wait for the job to finish, the results will be available in a zipped folder.

>Note: The free tier on `JuliaHub` offers 20 hours of computational time per month.

## Contributing

We welcome contributions in various forms, including bug reports, documentation improvements, feature suggestions, and more.

## Questions
For any questions or inquiries about PeriLab.jl, feel free to reach out to the authors via email.

## Authors and acknowledgment
[Dr.-Ing. Christian Willberg](mailto::christian.willberg@dlr.de)

[M.Sc. Jan-Timo Hesse](mailto::jan-timo.hesse@dlr.de)

## Project status
`PeriLab` is currently in development.
