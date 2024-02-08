# Getting Started

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

The simplest way to run the `PeriLab` simulation core is to use a provided example import the module and go. 

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
$ mpiexecjl -n 2 julia -e "using PeriLab; PeriLab.main()" Dogbone.yaml -v
```

## Training

The training input is given under the examples folder. The documentation and a video will follow.

## Index
```@index
Pages = ["gettin_started.md"]
```

## Functions
```@meta
CurrentModule = PeriLab
```

```@docs
main
get_examples
```
