<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# Changelog

All notable changes to this project will be documented in this file.

## [1.0.6] - 2024-

### Added

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

[1.0.2]: https://gitlab.com/dlr-perihub/PeriLab.jl/-/compare/v1.0.1...v1.0.2
[1.0.1]: https://gitlab.com/dlr-perihub/PeriLab.jl/-/compare/v1.0.0...v1.0.1
[1.0.0]: https://gitlab.com/dlr-perihub/PeriLab.jl/-/tags/v1.0.0
