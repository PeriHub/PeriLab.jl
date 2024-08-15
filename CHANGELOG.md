<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# Changelog

All notable changes to this project will be documented in this file.

## [1.2.1] - 2024-08-14

### Fixed

- Varying material properties #175

### Changed

- pre-commit hooks
- code cleaning

## [1.2.0] - 2024-08-09

### Added

- Bond associated formulation #154
- Surface extrusion for Abaqus
- Strain for non correspondence models

### Fixed

- Quasi-contact
- Inter critical value

## [1.1.6] - 2024-06-24

### Added

- Bond associated basis
- Many tests
- Element support for exodus export #20
- State variable in exodus export #148

### Fixed

- FEM basis, not yet coupled with PD #132
- Testset continues running after error

### Changed

- Moved datamanager #61
- Code cleaning and test coverage
- Optimized handling with specific volume

## [1.1.5] - 2024-05-13

### Fixed

- CSV Output order
- Force Boundary Condition
- Nodeset compute class
- Orientation
- UMAT Interface

### Added

- Solver summary
- Calculate cauchy and von Mises stress
- calculate_shape_tensor and calculate_deformation_gradient functions #152

### Changed

- Optimized Anisotropic Damage
- Optimized MPI communication #151

## [1.1.4] - 2024-04-10

### Fixed

- Exodus global export #145
- UMAT Header #139
- Anisotropic Damage #136

### Added

- UMAT #138
- Abaqus mesh to txt
- Docs

### Changed

- MPI non-blocking
- Bond length seperated from field #137

## [1.1.3] - 2024-03-28

### Fixed

- Volume calculation for Abaqus
- #137

### Added

- Abaqus Test

## [1.1.2] - 2024-03-25

### Fixed

- Abaqus surface extension

## [1.1.1] - 2024-03-21

### Fixed

- Abaqus input, still waiting for [Pull Request](https://github.com/JuliaFEM/AbaqusReader.jl/pull/71)

## [1.1.0] - 2024-03-19

### Changed

- Julia version

### Fixed

- MPI Issue
- Logging


## [1.0.7] - 2024-03-12

### Added

- Corrosion
- PrettyTables logging
- Logging Datetime
- MPI summary
- StaticArrays for performance

### Fixed

- Global export in MPI
- CSV export

### Changed

- Neighborhood Distribution
- CompactTension Example
- Docs

## [1.0.6] - 2024-02-14

### Fixed

- Silent mode allows log file

### Changed

- get_field uses an initialized function, performance improvements
- Allow older dependencies
- Dockerfile

## [1.0.5] - 2024-02-07

### Added

- Bond Filter Contact
- Surface Extension
- Von Mises Calculation
- Basic FEM support (#137, #136, #135, #134, #132, #129, #128, #124, #123, #122)
- Aqua test
- Variable datafield input (#139)
- Correspondence flexible material (#138)
- Plasticity (#120)
- Git info to logging

### Fixed

- Memory Leaks
- Reimport warnings

### Changed

- Read nodeset moved to core 1 (#140)

## [1.0.4] - 2024-01-11

### Added

- Abaqus mesh input (.inp)

## [1.0.3] - 2024-01-08

### Added

- Specific volume for additive models
- Tests for additive models
- Exodus input

### Fixed

- Memory issues with exodus

## [1.0.2] - 2023-12-18

### Added

- Anistropic Damage Model

## [1.0.1] - 2023-12-06

### Added

- JuliaHub support
- Codecov Support (#37)
- Docs
- PackageCompiler

### Fixed

- Wrong output path

### Changed

- Readme
- License
- Project.toml
- Optimized pd_solid
- Haskey with get()

### Removed

- Unecessary functions

## [1.0.0] - 2023-11-30

### Added

- First full PeriLab release

[1.0.2]: https://github.com/PeriHub/PeriLab.jl/-/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/PeriHub/PeriLab.jl/-/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/PeriHub/PeriLab.jl/-/tags/v1.0.0
