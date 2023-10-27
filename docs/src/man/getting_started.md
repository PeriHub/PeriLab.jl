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

## `PeriLab` 

The simplest way to run the 
`PeriLab` simulation core is to use a provided example import the module and go. 

```julia PeriLab
julia> using PeriLab

julia> PeriLab.main("examples/Dogbone/Dogbone.yaml")
```
The main functionalities for the `yaml` input deck is given in
```
"examples/functionalities.yaml"
```

```@meta
CurrentModule = PeriLab
```

```@docs
main
```
