# First Steps with PeriLab

## Setting up the Environment

PeriLab needs Julia 1.12 or newer. Install it as an app:
You can do it using the following commands:

```julia
using Pkg
Pkg.Apps.add("PeriLab")
```

or

```julia
julia> ] # ']' should be pressed

(@v1.12) pkg> app add PeriLab
```
`PeriLab --help` lists the command line options; the version is printed at the start of every run.

If you want to make sure everything works as expected you can run the tests
bundled with PeriLab:

```julia
using Pkg

Pkg.test("PeriLab") # This will take a few minutes.
```
