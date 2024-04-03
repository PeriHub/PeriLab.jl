## Development Guide

PeriLab is design to be extentable. Modules can be integrated with the so called factory modules. These modules are the interfaces to the higher functions. 

!!! info As an example: The solver calls compute_forces as a function. This function is integrated in the Material module. Within this module the relevant material model is integrated, using the [set_modules module](@ref "Module integration")