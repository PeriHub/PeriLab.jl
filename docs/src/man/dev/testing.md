# Software testing
The testing takes place manually or if commits are submitted to the repository. Make sure that all test are reasonably complex and won't take to much time. The testing [coverage](https://app.codecov.io/gh/perihub/PeriLab.jl) is checked as well.

!!! warning "Code quality"
    Tests only handle expected results and helps that code stays as intended after further development. It does not guarantee that not error occur.

!!! info "New tests"
    If new unexpected errors occur, please add an test to check possible exception to avoid them next time.

## Runtests.jl

The tests to perform are defined in the **runtests.jl** within the test folder. The structure is as follows.


```julia
@testset ExtendedTestSet "PeriLab" begin

    @testset "unit_tests" begin
        @testset "Compute" begin
            @testset "ut_compute_global_values" begin
                include(filepath + filename)
            end
            @testset "ut_compute_field_values" begin
                include(filepath + filename)
            end
        end
    end

    @testset "full_scale" begin
        ...
    end
```
After testing you will get an overview over failed or errored tests. The test set structure helps to identify the position of errors.


!!! info "Error & Fails"
    Error means that the code of the test has an error. Typically the function call is flawed. Fail means that the test is wrong, e.g.
    2==3



## Unit tests
Unit tests are tests to check functions. To do that you add a file the the unit_tests folder. The folder has the same structure as the src folder. Please add your test in place of the src code file where your function is located. It helps to find the tests. Add your tests (please take examples from already existing tests).


!!! info "Exceptions"
    In this tests the exeptions should be tested, because full scale tests are not capable of.


## Full scale test
Full scale models should be used to test features within a complete analyses. A complete model is run and tested against the result file.

**How to set up?**
To create such test, you have to create a ''normal'' model and run it. Create a folder of your test and put all the information for the model inside. Create .jl file in the following structure
```julia
    folder_name = basename(@__FILE__)[1:end-3]
    cd("fullscale_tests/" * folder_name) do
        run_perilab("additive_2d", 1, true, folder_name)
        run_perilab("additive_2d_heat", 1, true, folder_name)
        run_perilab("additive_3d", 1, true, folder_name)
    end
```

This file has to be called in the **runtests.jl**. You are able to define  .cmd with the following input if specific tolerances are needed.

```
DEFAULT TOLERANCE absolute 1.0E-9
COORDINATES absolute 1.0E-12
TIME STEPS absolute 1.0E-14
NODAL VARIABLES absolute 1.0E-12
    Temperature   absolute 1.0E-8
    Heat Flow     absolute 1.0E-8
    Active        absolute 1.0E-8
```


!!! warning "Naming convention"
    Names of the .yaml input and the output must be the same.


!!! info "License"
    Use licensing files for files which are not ascii.

## Examples
Examples should be given for more complex models used in papers. Here it should be refered to the commit hash (when it worked) and the paper where the model was used. These models won't run automatically. But they help to reproduce the results and make them more transparent.
