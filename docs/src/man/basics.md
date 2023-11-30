# First Steps with PeriLab

## Setting up the Environment

If want to use the PeriLab package you need to install it first.
You can do it using the following commands:

```julia
using Pkg

Pkg.add("PeriLab")
```

or

```julia
julia> ] # ']' should be pressed

(@v1.9) pkg> add PeriLab
```
Additionally, it is recommended to check the version of PeriLab that
you have installed with the `status` command.

```julia
julia> ]

(@v1.9) pkg> status PeriLab
      Status `~\v1.6\Project.toml`
  [a93c6f00] PeriLab v1.0.0
```

Throughout the rest of the tutorial we will assume that you have installed the
PeriLab package and have already typed `using PeriLab` which loads the
package:

```julia
using PeriLab
```
If you want to make sure everything works as expected you can run the tests
bundled with PeriLab, but be warned that it will take more than 30
minutes:

```julia
using Pkg

Pkg.test("PeriLab") # Warning! This will take more than 30 minutes.
```
