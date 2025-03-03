<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# Contributing
> If it is the first time that you contribute, please add yourself to the list
> of contributors below.

Contributions are always welcomed!

If you want to provide a code change, please:

* Create a clone or fork of the GitHub project.
* Develop the feature/patch.
* Provide a pull request.

> Take a look at the [GitHub Documentation](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) if you don't know how to do that.

## Contributors
* Christian Willberg
* Jan-Timo Hesse
* Felix Winkelmann
* Anna Pernatii
* Francisco A. Yapor Genao

## What can I contribute?

We are looking for different kind of contributions:

* Requirements (e.g. system requirements, installation instructions, etc.)
* Material models
* Solver implementations
* Performance improvements
* Documentation
* Testing (e.g. unit tests, integration tests, etc.)
* Mesh wrapping

## Code conventions

### Code style

* Make sure that [pre-commit](https://pre-commit.com/#installation) is installed and run `pre-commit install` before you commit.
* Functions and variables should be snake_case, e.g. `my_variable`.
* Functions and variables should clearly describe their purpose and function, e.g. `bond_forces` and not `b_f`.
* Modules should use capital letter `My_Module`
* Maximum of four if statement or loop level. If more are needed work with functions, return to exit, or functions like break or continue
* Type should be specified
* The integer and Float standard are Int64 and Float64
* Variables in headers of functions should have a type. If multiple types are possible use the Union{} option
* Input parameter check only in the initialization. These checks cost time during the solving process.

### Code quality

* Add tests for new features.
* Add your input values in the yaml check in parameter_handling.jl.
* Do not call Modules at same level, e.g. Material Models should not call Thermal Models.

### Documentation

* Add documentation for new features.

## Code Review

The following reviewers will review your pull request. The duration may vary depending on the size of the pull request.

We are trying to keep the review process as short as possible.

### Reviewers
* Christian Willberg
* Jan-Timo Hesse
