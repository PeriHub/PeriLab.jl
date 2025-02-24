# Parameters
All the parameter are stored in a dictionary called params.

The essential parameter are checked in [parameter_handling.jl](https://github.com/PeriHub/PeriLab.jl/src/support/parameters/parameter_handling.jl) with the variable global _expected_structure_. If parameter are not defined there, a warning is given in the log file.

Parameter are included via the YAML input deck. The parameter can be found in the program in a dictionary. Therefore, the parameter naming is used there as well. For some variables, e.g. boundary condition, equations can be specified. These equations are interpreted and can be used to include time or spatial depended variables.

!!! note "Good start"
    Please check some of the full scale tests. There are several yaml files with parameter definitions.
