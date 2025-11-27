<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# Changelog

All notable changes to this project will be documented in this file.

## [1.5.1] - 2025-11-27

### Added

- Exodus SPHERE element input type

### Fixed

- Cancellation for max. damage

## [1.5.0] - 2025-10-30

### Added

- PeriLabError, to handle known exceptions
- FieldTypes added and field creation streamlined

### Changed

- Data_Manager removed as an argument
- Removed module includes, reducing compilation time
- Module names streamlined
- TimerOutputs is now using the default timer


## [1.4.10] - 2025-10-06

### Added

- Result files will be merged even after an error ocurs

### Changed

- Reduced allocations via new function interface
- Optimized pre calculation factory
- Docker image updated

## [1.4.9] - 2025-08-13

### Added

- Print Bed Z Coordinate
- Time dependent environmental temperature
- Allow state variables for compute parameters

### Fixed

- Optimized Print Bed thermal flow

## [1.4.8] - 2025-08-07

### Added

- Curved Gcode interpretation
- Nearest Point Data Compute Class
- Global Variable fullscale test
- Fieldtype check befor export #287

### Fixed

- Global Variable calculation for MPI #285
- Nodeset evaluation for 2D

### Changed

- Moved the anisotropic damage calculation to a new material model

## [1.4.7] - 2025-07-17

### Fixed

- Bond rotation for heat transfer

## [1.4.6] - 2025-07-15

### Added

- Dwell (G4) callback #292
- 1D Bond based example

### Fixed

- Bond field export #288
- Node set id assignment
- Global value export for MPI

### Changed

- Reworked the gcode reader and specifc volume
- Optimized heat transfer module
- Substructured the datamanager for better readability
- Improve typstability and MPI calls

## [1.4.5] - 2025-06-30

### Added

- Only surface nodes options for contact
- Contact search frequency #279
- Documentation for search strategy #282

### Fixed

- Contact penalty model
- Contact for MPI #280
- Type stability

### Changed

- Contact search algorithm #281

## [1.4.4] - 2025-06-23

### Added

- Start, Stop and End command definiton for gcode input

## [1.4.3] - 2025-06-20

### Fixed

- Fixed parameter handling of globals #276
- Fixed a bug with the gcode import

## [1.4.2] - 2025-06-19

### Added

- Friction Coefficient for contact model #275
- Damage initiation flag
- Error Handling for eval of node sets
- Contact groups for multiple use of the same contact model

### Fixed

- Mesh angle definition
- MPI for conatct models

### Changed

- Contact search optimized
- BC evaluation optimized with function call
- Fields for all steps are initialized at first step
- Optimized plastic material model

## [1.4.1] - 2025-05-21

### Added

- Contact search radius

### Fixed

- Logging for silent mode
- Consisten Surface IDs #268

## [1.4.0] - 2025-05-14

### Added

- Multibody contact support #264
- cos, sin, tan BC functions for vector operations

### Fixed

- Contact Surface Detection
- Rotation tensor switch
- Logging path for filenames with underscores in the name

### Changed

- Improved critical stretch performance

## [1.3.7] - 2025-04-22

### Added

- Thermal decomposition model
- Dependent values for inter block damage model
- Block summary now also includs if a block is PD or FEM

### Fixed

- Block ID definition for MPI
- Static solver for non-thermal analysis
- BC check for missing variables
- STATEV exodus export for UMAT and VUMAT

### Changed

- Renamed corrosion to degredation
- Number of steps for static solver
- Nodeset definition, now supports x, y and z

## [1.3.6] - 2025-04-01

### Added

- [SciML Style](https://github.com/SciML/SciMLStyle)
- Reference temperature for the thermal expansion model #252
- Thermal expasion is now also usable with the static solver #248

### Fixed

- Static solver for multistep simulation

### Changed

- Neighborhood search is merged into the 'get_nearest_neighbors' function

## [1.3.5] - 2025-03-19

### Added

- Additional Time parameter for multistep analysis
- Block and nodeset definition via input deck
- Gcode files as input mesh

### Fixed

- Bond force calculation for correspondence models
- Init fields will only be called once

## [1.3.4] - 2025-03-10

### Added

- Surface correction #240
- Static solver start solution
- Static to verlet switch #238

### Fixed

- Critical stretch damage model
- Static solver
- Global vlues clculation
- VUMAT without temperature models

### Changed

- Boundary cnditions function

## [1.3.3] - 2025-02-18

### Added

- Multistep Solver #231
- VUMAT interface
- Documentations #232 #229 #228

### Fixed

- Critical time step estimation for orthotropic materials
- "Safety factor" warning

### Changed

- Deformation gradient is only calculated when needed
- PeriLab Banner width

## [1.3.2] - 2025-02-05

### Fixed

- Hooke matrix calculation for orthotropic materials
- Zero Energy Control with rotations #226
- HEATVAL interface speed

## [1.3.1] - 2025-02-02

### Breaking Changes

- Added "Block ID" definition

### Added

- Start and End Time for Output definition

### Changed

- Block name can be anything now
- Abaqus element block definition
- Renamed to "Calculate von Mises stress"

### Fixed

- UMAT interface parameter definition and speed
- Correspondence Plastic #223

## [1.3.0] - 2025-01-24

### Added

- FEM-PD Coupling #218
- Multiple Job Excution #222

## [1.2.7] - 2025-01.20

### Added

- Temperature dependent orthotropic material #215
- Blockwise angle definition

### Fixed

- Field allocation for additive
- MPI Testing

## [1.2.6] - 2025-01-10

### Added

- Unified Bond-based Elastic model #210
- Orthothropic material definition for Correspondence models #216
- Gcode reader script
- Dependency properties for Energy Release Damage Model and Bond-based Elastic #215

### Changed

- Bond fields are saved as vectors, reducing memory usage
- Optimzed memory allocation and speed
- Optimized additve process simuations

## [1.2.5] - 2024-10-24

### Fixed

- HETVAL subroutine
- Correspondence routines
- Damage index function

### Changed

- Optimzed memory allocation and speed
- Refine documentations

## [1.2.4] - 2024-10-08

### Added

- Include lambda rotation in Thermal flow #134
- Int support for solver params

### Fixed

- Many performance improvements
- Fix high abaqus mesh memory usage
- Fix OutOfMemory MPI Error

### Changed

- Removed view from get_field method and streamlined it

## [1.2.3] - 2024-09-13

### Added

- Calculation test

### Fixed

- Calculation methods

## [1.2.2] - 2024-09-13

### Added

- HEATVAL fortran routine
- Two block test #195

### Fixed

- Test cmd config files
- calculate_nodelist
- Logging
- Fix Poisson's ratio #196

### Changed

- Physics module is now the Model_Factory #198
- Inputdeck: Definition changes from Physics to Models
- Restructuring of the physics factory #168
- Block Wise Defintion of pre calculation #193
- Move pre calculation in models #191
- Optimized damage models and geometry routines #188

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

- Degradation
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
