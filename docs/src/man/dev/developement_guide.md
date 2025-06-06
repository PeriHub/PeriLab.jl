## Development Guide

PeriLab is design to be extentable. Modules can be integrated with the so called factory modules. These modules are the interfaces to the higher functions.

!!! info "Example"
    The solver calls compute\_forces as a function. This function is integrated in the Material module. Within this module the relevant material model is integrated, using the [set_modules module](@ref "Module integration")

## Model structure

PeriLab is structured in different ''Models''. These models are included in the model factory and split in two main parts. An init and compute part. All memory, parameter check, etc. takes place in the init phase. The computation should run without thowing warnings or errors due to input variable errors.

If you plan large scale implemention please contact us.

## Guidelines
- try to fix the type of the variables in functions
- use self explaining variables and functions
- write tests
    - full scale tests should check whole models and their interaction in the software
    - unit tests should test exceptions and key features
- inputs are checked in the init part of models not during run time
- in code documentation in front of functions
- theory and extensive model explainations in the documentation
- write issues
- ask for help
- have fun :)

!!! info "References"
    If you are part of the developer or user of PeriLab feel free to use your references to improve the documentation. Especially for publications created with PeriLab.
