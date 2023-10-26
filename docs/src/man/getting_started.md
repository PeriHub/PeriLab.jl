# Getting Started

## Installation

The DataFrames package is available through the Julia package system and can be installed using the following commands:

```julia
using Pkg
Pkg.add("PeriLab")
```

Throughout the rest of this tutorial, we will assume that you have installed the
PeriLab package and have already typed `using PeriLab` to bring all of the
relevant variables into your current namespace.

## The `DataFrame` Type

Objects of the `DataFrame` type represent a data table as a series of vectors,
each corresponding to a column or variable. The simplest way of constructing a
`DataFrame` is to pass column vectors using keyword arguments or pairs:

```julia PeriLab
julia> using PeriLab

julia> PeriLab.main("examples/Dogbone/Dogbone.yaml")
```

```@meta
CurrentModule = PeriLab
```

```@docs
main(filename, dry_run, verbose, debug, silent)
PeriLab
```
