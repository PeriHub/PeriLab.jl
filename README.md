<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# `PeriLab` - Peridynamic Laboratory

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://perihub.github.io/PeriLab.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://perihub.github.io/PeriLab.jl/stable

[ci-img]: https://github.com/perihub/PeriLab.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/perihub/PeriLab.jl/actions/workflows/CI.yml

[cov-img]: https://codecov.io/gh/perihub/PeriLab.jl/branch/main/graph/badge.svg
[cov-url]: https://codecov.io/gh/perihub/PeriLab.jl

[code-style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[code-style-url]: https://github.com/invenia/BlueStyle

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[sciml-img]: https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826
[sciml-url]: https://github.com/SciML/SciMLStyle

[doi-img]: https://zenodo.org/badge/DOI/10.1016/j.softx.2024.101700.svg
[doi-url]: https://doi.org/10.1016/j.softx.2024.101700

[release-img]: https://img.shields.io/github/v/release/PeriHub/PeriLab.jl
[release-url]: https://github.com/PeriHub/PeriLab.jl/releases

[docker-img]: https://img.shields.io/docker/pulls/perihub/perilab
[docker-url]: https://hub.docker.com/r/perihub/perilab

[license-img]: https://img.shields.io/badge/License-BSD-blue.svg
[license-url]: https://github.com/PeriHub/PeriLab.jl/LICENSE

[youtube-img]: https://img.shields.io/youtube/channel/subscribers/UCeky7HtUGlOJ2OKknvl6YnQ
[youtube-url]: https://www.youtube.com/@PeriHub

[hirse-img]: https://img.shields.io/badge/Promo-8db427?label=HiRSE&labelColor=005aa0&link=https%3A%2F%2Fgo.fzj.de%2FCodePromo

| **Documentation** | **Build Status** |  **Quality** |
|:----:|:----:|:----:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ci-img]][ci-url] [![][cov-img]][cov-url] | [![][aqua-img]][aqua-url] [![][sciml-img]][sciml-url] |
| **Deployment** | **License** | **Socials** |
| [![][release-img]][release-url]  [![][docker-img]][docker-url] | [![][license-img]][license-url] [![][doi-img]][doi-url] | [![][youtube-img]][youtube-url] ![HiRSE Code Promo 
Badge][hirse-img] |

Welcome to `PeriLab`, a comprehensive platform for working with peridynamics, designed to support researchers and practitioners at all levels of expertise.

<p align="center" style="font-size:0;"><!--
  PeriLab_crack      --><img align="middle" src="https://raw.githubusercontent.com/PeriHub/PeriLab.jl/main/assets/PeriLab_crack.gif" width="50%"><!--
  PeriLab_additive      --><img align="middle" src="https://raw.githubusercontent.com/PeriHub/PeriLab.jl/main/assets/PeriLab_additive.gif" width="50%">
</p>

## Overview
- [Features](#features)
- [Documentation](#documentation)
   - [Examples](#examples)
- [Installation](#installation)
   - [Getting Started with `PeriLab`](#getting-started-with-perilab)
   - [Parallel Processing with PeriLab (MPI)](#parallel-processing-with-perilab-mpi)
   - [Installing with Docker](#installing-with-docker)
   - [Web Framework PeriHub](#web-framework-perihub)
   - [PeriLab on PeriHub](#perilab-on-perihub)
- [What's Next](#whats-next)
- [Contributing](#contributing)
   - [Development](#development)
- [Questions](#questions)
- [How to cite](#how-to-cite)
- [Acknowledgments](#acknowledgments)
   - [Partners](#partners)
   - [Authors](#authors)
- [Project status](#project-status)
## Features

- 🚀 **Easy Installation**: PeriLab's straightforward installation process makes it accessible for researchers and engineers without extensive computational expertise.

- ✒️ **Modularization**: The software is designed with a modular architecture that allows users to easily integrate their own material and damage models.

- 🎨 **Formulations**: Bond-based, bond-associated, as well as oridnary and non-ordinary state-based peridynamic formulations can be used with PeriLab.

- 🔩 **Material models**: PeriLab supports various material models, such as elastic, plastic, and more, enabling simulation of complex materials and structures.

- 🔨 **Damage models**: [Damage](https://www.youtube.com/watch?v=ClV2ojQPrFM) models such as critical stretch or an energy based criterium are included to simulate different types of damage, such as crack propagation or delamination, in their peridynamic simulations.

- 👊 **Contact**: [Simulation](https://www.youtube.com/watch?v=qj7xGgmjEdE) of contact between objects is supported with a variety of contact conditions.

- 👬 **FEM/PD coupling**: Coupling between Peridynamics and the Finite element method is supported.

- 🔥 **Additive Manufacturing**: PeriLab supports [additive manufacturing](https://www.youtube.com/watch?v=J55n2xWJosA), allowing users to create custom additive models for their simulations.

- 🔑 **Solver**: Different solvers like `Verlet` or `Static` can be utilized, swiching between those is also supported.

- 🧲 **Multi physics**: PeriLab supports [multimodels](https://www.youtube.com/watch?v=byhz-_RW0Lw) simulations, combining different types of peridynamic physical  and damage models to create a comprehensive simulation environment.

- 🔌 **Subroutine Interfaces**: Subroutines, such as UMAT, VUMAT and HETVAL are usable as material models

- ⚡ **MPI**: PeriLab supports parallel computing using Message Passing Interface (MPI) technology to improve simulation performance on high-performance clusters.

- 🔁 **Multistep simulations**: PeriLab supports the definition of multiple solver steps, allowing to combine different environmental conditions in a single run.

- 📐 **Surface Correction**: PeriLab provides tools for surface correction, such as the `Volume correction` method.

- 💻 **HPC capabilities**: PeriLab is designed for high-performance computing (HPC) environments, allowing users to run large-scale simulations efficiently.

- 📤📥 **Exodus Input/Output**: PeriLab uses the Exodus II data format for input and output, enabling easy transfer of data between simulation tools.

- 🧮 **Abaqus Input**: PeriLab supports Abaqus input files, allowing users to create custom Abaqus models for their simulations.

- ➗ **Bond filter**: The bond filter feature allows users to apply specific conditions to the bonds between particles in a simulation, influencing their behavior and interaction with other particles.

- 🔧 **User specified Input/Output**: PeriLab provides flexible options for users to specify custom input and output files, enabling easy data manipulation and analysis.

- 🧪 **Test Pipeline**: The PeriLab Source Code will be tested in a test pipeline to ensure its correctness and performance.


## Documentation

Explore the comprehensive [documentation](https://perihub.github.io/PeriLab.jl/) for `PeriLab` and
the seminar information for the first German Peridynamics course [Lecture Non-local structural mechanics and peridynamics](https://perihub.github.io/PeriLab.jl/dev/lecture/lecture/).

### Examples

A few basic examples of `PeriLab` can be found in the [examples](https://github.com/PeriHub/PeriLab.jl/tree/main/examples) directory, or if you want to have a look at results go to our growing [PeriLab-Results service](https://perilab-results.nimbus-extern.dlr.de).


## Installation

The `PeriLab`  package is available through the Julia package system and can be installed using the following commands:

```julia
using Pkg
Pkg.add("PeriLab")
```

Throughout the rest of this tutorial, we will assume that you have installed the
PeriLab package and have already typed `using PeriLab` to bring all of the
relevant variables into your current namespace.

### Getting Started with `PeriLab`

Jumpstart your exploration of the PeriLab simulation core with provided examples. Run the following commands in Julia:

```julia PeriLab
using PeriLab

PeriLab.get_examples()
PeriLab.main("examples/DCB/DCBmodel.yaml")
```
>Note: More details about the main functionalities in the yaml input deck [here](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Support/Parameters/parameter_handling.jl).

### Parallel Processing with `PeriLab` (MPI)

To handle large-scale problems efficiently, install [MPI](https://juliaparallel.org/MPI.jl/stable/usage/). Run PeriLab with two processors on a **Linux** system using the following commands:

```sh
$ julia
julia> using MPI
julia> MPI.install_mpiexecjl()
```
>Note: If you work with **Windows 10 or higher** you can use the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) environment.

Run PeriLab with two processors:
```sh
$ mpiexecjl -n 2 julia -e 'using PeriLab; PeriLab.main("examples/DCB/DCBmodel.yaml")'
```

>Note: For HPC configurations please refer to [here](https://juliaparallel.org/MPI.jl/stable/configuration/#configure_jll_binarys).

### Installing with Docker

 To install PeriLab using the official Perihub/Perilab Docker image, follow these steps:

1. **Install Docker**: Before you begin, ensure that you have Docker installed on your system. You can download and install Docker from the official website (https://www.docker.com/). Make sure your system meets the minimum requirements for running Docker.

2. **Pull the Perihub/Perilab Docker image**: Use the following command in a terminal or command prompt to pull the latest Perihub/Perilab Docker image from the Docker Hub repository:

   ```bash
   docker pull perihub/perilab
   ```

3. **Run the Docker container**: Once the image has been downloaded, create a new directory for your PeriLab simulations and navigate to it in the terminal or command prompt. Run the following command to start the Docker container:

   ```bash
   docker run -it --rm -v <path_to_local_simulations_directory>:/app/simulations perihub/perilab bash
   ```

   Replace `<path_to_local_simulations_directory>` with the absolute path to a local directory where you want to store your PeriLab simulations. This command will open a new terminal session inside the Docker container.

Now, you've successfully installed PeriLab using the official Perihub/Perilab Docker image. You can start running your own peridynamic simulations within the container.

### Web Framework `PeriHub`

PeriLab is also included as a ready to use application in the [PeriHub](https://github.com/PeriHub/PeriHub) web framework.

### `PeriLab` on `JuliaHub`

Experience the convenience of using PeriLab as a ready-to-use application on JuliaHub. Simply create an [account](https://juliahub.com), navigate to the [applications page](https://juliahub.com/ui/Applications), and add the repository URL: https://github.com/PeriHub/PeriLab.jl.

Configure advanced options, such as _filename_, _dryrun_, _verbosity_, _debug_, and _silence_. Click __Start__ and monitor the job progress. Results will be available in a zipped folder.

Hit the __Start__ button and wait for the job to finish, the results will be available in a zipped folder.

>Note: The free tier on `JuliaHub` offers 20 hours of computational time per month.

## What's Next?

Here are some exciting tasks on our roadmap:

- 💧 **Degradation**: The simulation of degrative materials will be added to PeriLab

- ➕ **More material and damage models**: PeriLab's future development plans include adding more sophisticated material models (e.g., viscoelastic-plastic) and damage models, expanding the software's applicability to a wider range of real-world phenomena.

- ✂️ **Distribution logic**: As part of its ongoing development, PeriLab will continue to incorporate new distribution logic for improved performance and reduced computational resources.

- 🏎️ **Optimizations**: As part of its ongoing development, PeriLab will continue to focus on optimizing the simulation process by incorporating new techniques like parallel optimization algorithms for improved efficiency and reduced computational resources.

Feel free to contribute and help us make PeriLab even better! 🙌

Here are the planned [issues](https://github.com/PeriHub/PeriLab.jl/issues) and [milestones](https://github.com/PeriHub/PeriLab.jl/milestones).

## Contributing

We welcome contributions in various forms, including bug reports, documentation improvements, feature suggestions, and more. To get started, follow these steps and have a look at the [Contribution Guidelines](CONTRIBUTING.md):

### Development
1. **Clone the repository:**
```sh
git clone https://github.com/PeriHub/PeriLab.jl
cd PeriLab.jl
```
2. **Activate the environment and install dependencies:**
```sh
$ julia
julia> ]
pkg> activate .
pkg> up
```
3. **Run the script:**
```sh
$ julia --project=. src/main.jl examples/DCB/DCBmodel.yaml
```

## Questions
For any questions or inquiries about PeriLab.jl, feel free to reach out to the authors via E-mail or open a [discussion](https://github.com/PeriHub/PeriLab.jl/discussions).


## How to cite
To cite PeriLab in your publications please use the following [paper](https://doi.org/10.1016/j.softx.2024.101700).

```s
@Article{WillbergC2024,
author={Willberg, Christian
and Hesse, Jan-Timo
and Pernatii, Anna},
title={{PeriLab - Peridynamic Laboratory}},
journal={SoftwareX},
year={2024},
publisher={Elsevier},
volume={26},
issn={2352-7110},
doi={10.1016/j.softx.2024.101700},
url={https://doi.org/10.1016/j.softx.2024.101700}
}
```





## Acknowledgments
Funding and acknowledgment infos can be found under [acknowledgments](ACKNOWLEDGEMENTS.md)

### Partners


| <img src="https://raw.githubusercontent.com/PeriHub/PeriLab.jl/main/assets/dlr.jpg" height="200" title="German Aerospace Center"> |
:------------------------------------------------------------------------------------------------------------------------------:|
| [German Aerospace Center](http://www.dlr.de/sy)                                                                               |
|<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/5/5c/Logo_h2.svg/1280px-Logo_h2.svg.png" height="200" title="Magdeburg-Stendal University of Applied Science"> |
[Magdeburg-Stendal University of Applied Science](http://www.h2.de)                                                                  |
 <img src="https://www.cd.ovgu.de/cd_media/CD_OVGU/Downloads/Logo_jpg_png_svg_EPS_pdf/Logodownload/OVGU_Logo-download-1-p-1960.jpeg" height="200" title="Otto von Guericke University Magdeburg"> |
|[Otto von Guericke University Magdeburg](http://www.ovgu.de)|
 <img src="https://cdn.chimpify.net/5b34fb88a8587269268b457d/2018/07/SWMS_logo.png" height="200" title="SWMS Systemtechnik Ingenieurgesellschaft mbH"> |
|[SWMS Systemtechnik Ingenieurgesellschaft mbH](https://www.swms.de)|

### Authors
<p>

<a href="https://orcid.org/0000-0003-2433-9183"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" style="height:15px;width:auto;vertical-align: top;background-color:transparent;"> </a>[Prof. Dr.-Ing. Christian Willberg](mailto::christian.willberg@h2.de)

</p>

<p>

<a href="https://orcid.org/0000-0002-3006-1520"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg"  style="height:15px;width:auto;vertical-align: top;background-color:transparent;"> [M.Sc. Jan-Timo Hesse](mailto::jan-timo.hesse@dlr.de)

</p>

## Project status
`PeriLab` is currently in development.
