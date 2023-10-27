# PeriLab

Welcome to the PeriLab documentation!

This resource aims to teach you everything you need to know to get up and
running with tabular data manipulation using the PeriLab.jl package.

For more illustrations of PeriLab.jl usage, in particular in conjunction with
other packages you can check-out the following resources
(they are kept up to date with the released version of PeriLab.jl):
* [JuliaCon 2019](https://github.com/bkamins/JuliaCon2019-PeriLab-Tutorial),
 
## What is PeriLab.jl?

PeriLab.jl provides a Peridynamics [BobaruF2016](@cite) simulation core

Its design and functionality are similar to those of Peridigm  [LittlewoodDJ2023](@cite) and several extenions [WillbergC2019](@cite), [WillbergC2023](@cite), [HesseJA2023](@cite). 


## PeriLab.jl and the Julia Data Ecosystem

The Julia data ecosystem can be a difficult space for new users to navigate, in
part because the Julia ecosystem tends to distribute functionality across
different libraries more than some other languages. Because many people coming
to PeriLab.jl are just starting to explore the Julia data ecosystem, below is
a list of well-supported libraries that provide different data science tools,
along with a few notes about what makes each library special, and how well
integrated they are with PeriLab.jl.


- **Statistics**
    - [StatsKit.jl](https://github.com/JuliaStats/StatsKit.jl): A convenience
      meta-package which loads a set of essential packages for statistics,
      including those mentioned below in this section and PeriLab.jl itself.

While not all of these libraries are tightly integrated with PeriLab.jl,
because `DataFrame`s are essentially collections of aligned Julia vectors, so it
is easy to (a) pull out a vector for use with a non-PeriLab-integrated
library, or (b) convert your table into a homogeneously-typed matrix using the
`Matrix` constructor or StatsModels.jl.

## Questions?

If there is something you expect PeriLab to be capable of, but
cannot figure out how to do, please reach out with questions in Domains/Data on
[Discourse](https://discourse.julialang.org/new-topic?title=[PeriLab%20Question]:%20&body=%23%20Question:%0A%0A%23%20Dataset%20(if%20applicable):%0A%0A%23%20Minimal%20Working%20Example%20(if%20applicable):%0A&category=Domains/Data&tags=question).
Additionally you might want to listen to an introduction to PeriLab.jl on
[JuliaAcademy](https://juliaacademy.com/p/introduction-to-PeriLab-jl).

Please report bugs by
[opening an issue](https://github.com/JuliaData/PeriLab.jl/issues/new).

You can follow the **source** links throughout the documentation to jump right
to the source files on GitHub to make pull requests for improving the
documentation and function capabilities.

Please review [PeriLab contributing
guidelines](https://github.com/JuliaData/PeriLab.jl/blob/main/CONTRIBUTING.md)
before submitting your first PR!

Information on specific versions can be found on the [Release
page](https://github.com/JuliaData/PeriLab.jl/releases).

## Package Manual

```@contents
Pages = ["man/basics.md",
         "man/getting_started.md"]
Depth = 2
```

## API

Only exported (i.e. available for use without `PeriLab.` qualifier after
loading the PeriLab.jl package with `using PeriLab`) types and functions
are considered a part of the public API of the PeriLab.jl package. In general
all such objects are documented in this manual (in case some documentation is
missing please kindly report an issue
[here](https://github.com/JuliaData/PeriLab.jl/issues/new)).

!!! note

    Breaking changes to public and documented API are avoided in
    PeriLab.jl where possible.
