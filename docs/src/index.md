# PeriLab

Welcome to the PeriLab documentation!

This resource aims to teach you everything you need to know to get up and
running with tabular data manipulation using the DataFrames.jl package.

For more illustrations of DataFrames.jl usage, in particular in conjunction with
other packages you can check-out the following resources
(they are kept up to date with the released version of DataFrames.jl):
* [JuliaCon 2019](https://github.com/bkamins/JuliaCon2019-DataFrames-Tutorial),
  [JuliaCon 2020](https://github.com/bkamins/JuliaCon2020-DataFrames-Tutorial),
  [JuliaCon 2021](https://github.com/bkamins/JuliaCon2021-DataFrames-Tutorial),
  [JuliaCon 2022](https://github.com/bkamins/JuliaCon2022-DataFrames-Tutorial),
  [PyData Global 2020](https://github.com/bkamins/PyDataGlobal2020),
  and [ODSC Europe 2021](https://github.com/bkamins/ODSC-EUROPE-2021) tutorials
* [DataFrames.jl showcase](https://github.com/bkamins/DataFrames-Showcase)

If you prefer to learn DataFrames.jl from a book you can consider reading:
* [Julia for Data Analysis](https://github.com/bkamins/JuliaForDataAnalysis);
* [Julia Data Science](https://juliadatascience.io/).

## What is DataFrames.jl?

DataFrames.jl provides a set of tools for working with tabular data in Julia.
Its design and functionality are similar to those of
[pandas](https://pandas.pydata.org/) (in Python) and `data.frame`,
[`data.table`](https://rdatatable.gitlab.io/data.table/)
and [dplyr](https://dplyr.tidyverse.org/) (in R),
making it  a great general purpose data science tool

## DataFrames.jl and the Julia Data Ecosystem

The Julia data ecosystem can be a difficult space for new users to navigate, in
part because the Julia ecosystem tends to distribute functionality across
different libraries more than some other languages. Because many people coming
to DataFrames.jl are just starting to explore the Julia data ecosystem, below is
a list of well-supported libraries that provide different data science tools,
along with a few notes about what makes each library special, and how well
integrated they are with DataFrames.jl.


- **Statistics**
    - [StatsKit.jl](https://github.com/JuliaStats/StatsKit.jl): A convenience
      meta-package which loads a set of essential packages for statistics,
      including those mentioned below in this section and DataFrames.jl itself.

While not all of these libraries are tightly integrated with DataFrames.jl,
because `DataFrame`s are essentially collections of aligned Julia vectors, so it
is easy to (a) pull out a vector for use with a non-DataFrames-integrated
library, or (b) convert your table into a homogeneously-typed matrix using the
`Matrix` constructor or StatsModels.jl.

## Questions?

If there is something you expect DataFrames to be capable of, but
cannot figure out how to do, please reach out with questions in Domains/Data on
[Discourse](https://discourse.julialang.org/new-topic?title=[DataFrames%20Question]:%20&body=%23%20Question:%0A%0A%23%20Dataset%20(if%20applicable):%0A%0A%23%20Minimal%20Working%20Example%20(if%20applicable):%0A&category=Domains/Data&tags=question).
Additionally you might want to listen to an introduction to DataFrames.jl on
[JuliaAcademy](https://juliaacademy.com/p/introduction-to-dataframes-jl).

Please report bugs by
[opening an issue](https://github.com/JuliaData/DataFrames.jl/issues/new).

You can follow the **source** links throughout the documentation to jump right
to the source files on GitHub to make pull requests for improving the
documentation and function capabilities.

Please review [DataFrames contributing
guidelines](https://github.com/JuliaData/DataFrames.jl/blob/main/CONTRIBUTING.md)
before submitting your first PR!

Information on specific versions can be found on the [Release
page](https://github.com/JuliaData/DataFrames.jl/releases).

## Package Manual

```@contents
Pages = ["man/basics.md",
         "man/getting_started.md"]
Depth = 2
```

## API

Only exported (i.e. available for use without `DataFrames.` qualifier after
loading the DataFrames.jl package with `using DataFrames`) types and functions
are considered a part of the public API of the DataFrames.jl package. In general
all such objects are documented in this manual (in case some documentation is
missing please kindly report an issue
[here](https://github.com/JuliaData/DataFrames.jl/issues/new)).

!!! note

    Breaking changes to public and documented API are avoided in
    DataFrames.jl where possible.
