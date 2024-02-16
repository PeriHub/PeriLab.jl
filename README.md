<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# `PeriLab` - Peridynamic Laboratory
Welcome to `PeriLab`, a powerful software solution designed for tackling Peridynamic problems. 

<p align="center" style="font-size:0;"><!--
  PeriLab_crack      --><img align="middle" src="https://gitlab.com/dlr-perihub/PeriLab.jl/-/raw/main/assets/PeriLab_crack.gif" width="50%"><!--
  PeriLab_additive      --><img align="middle" src="https://gitlab.com/dlr-perihub/PeriLab.jl/-/raw/main/assets/PeriLab_additive.gif" width="50%">
</p>

## Documentation

Explore the comprehensive [documentation](https://dlr-perihub.gitlab.io/PeriLab.jl/) for `PeriLab`

## Features ⭐

- 🚀 **Easy Installation**: PeriLab's straightforward installation process makes it accessible for researchers and engineers without extensive computational expertise.

- ✒️ **Modularization**: The software is designed with a modular architecture that allows users to easily integrate their own material and damage models.

- 🔩 **Material models**: PeriLab supports various material models, such as elastic, plastic, and more, enabling simulation of complex materials and structures.

- 🔨 **Damage models**: Damage models such as critical stretch or an energy based criterium are included to simulate different types of damage, such as crack propagation or delamination, in their peridynamic simulations.

- 🔥 **Additive Manufacturing**: PeriLab supports additive manufacturing, allowing users to create custom additive models for their simulations.

- ⚡ **MPI**: PeriLab supports parallel computing using Message Passing Interface (MPI) technology to improve simulation performance on high-performance clusters.

- 💻 **HPC capabilities**: PeriLab is designed for high-performance computing (HPC) environments, allowing users to run large-scale simulations efficiently.

- 📤📥 **Exodus Input/Output**: PeriLab uses the Exodus II data format for input and output, enabling easy transfer of data between simulation tools.

- 🧮 **Abaqus Input**: PeriLab supports Abaqus input files, allowing users to create custom Abaqus models for their simulations.

- ➗ **Bond filter**: The bond filter feature allows users to apply specific conditions to the bonds between particles in a simulation, influencing their behavior and interaction with other particles.

- 🔧 **User specified Input/Output**: PeriLab provides flexible options for users to specify custom input and output files, enabling easy data manipulation and analysis.

- 🧪 **Test Pipeline**: The PeriLab Source Code will be tested in a test pipeline to ensure its correctness and performance.

## Installation

The `PeriLab`  package is available through the Julia package system and can be installed using the following commands:

```julia
using Pkg
Pkg.add("PeriLab")
```

Throughout the rest of this tutorial, we will assume that you have installed the
PeriLab package and have already typed `using PeriLab` to bring all of the
relevant variables into your current namespace.

## Getting Started with `PeriLab` 

Jumpstart your exploration of the PeriLab simulation core with provided examples. Run the following commands in Julia:

```julia PeriLab
using PeriLab

PeriLab.get_examples()
PeriLab.main("examples/DCB/DCBmodel.yaml")
```
>Note: More details about the main functionalities in the yaml input deck [here](https://gitlab.com/dlr-perihub/PeriLab.jl/-/blob/main/src/Support/Parameters/parameter_handling.jl?ref_type=heads).

## Parallel Processing with `PeriLab` (MPI)

To handle large-scale problems efficiently, install [MPI](https://juliaparallel.org/MPI.jl/stable/usage/). Run PeriLab with two processors using the following commands:

```sh
$ julia
julia> using MPI
julia> MPI.install_mpiexecjl()
```

Run PeriLab with two processors:
```sh
$ mpiexecjl -n 2 julia -e "using PeriLab; PeriLab.main()" examples/DCB/DCBmodel.yaml -v
```

## Installing with Docker 🐳

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

## `PeriLab` on `JuliaHub`

Experience the convenience of using PeriLab as a ready-to-use application on JuliaHub. Simply create an [account](https://juliahub.com), navigate to the [applications page](https://juliahub.com/ui/Applications), and add the repository URL: https://gitlab.com/dlr-perihub/PeriLab.jl.

Configure advanced options, such as _filename_, _dryrun_, _verbosity_, _debug_, and _silence_. Click __Start__ and monitor the job progress. Results will be available in a zipped folder.

Hit the __Start__ button and wait for the job to finish, the results will be available in a zipped folder.

>Note: The free tier on `JuliaHub` offers 20 hours of computational time per month.

## What's Next? 🚀

Here are some exciting tasks on our roadmap:

- 🔑 **Quasi-static solver**: A future development for PeriLab is extending its capabilities with a more robust quasi-static solver for larger systems and complex boundary conditions.

- 👊 **Contact**: An upcoming feature in PeriLab is enhancing contact modeling to support advanced features like friction, adhesion, and contact forces based on temperature or other variables.

- ➕ **More material and damage models**: PeriLab's future development plans include adding more sophisticated material models (e.g., viscoelastic-plastic) and damage models, expanding the software's applicability to a wider range of real-world phenomena.

- 👬 **FEM/PD coupling**: A future enhancement for PeriLab is improving its FEM/PD coupling functionality by implementing more advanced techniques, such as a seamless data exchange between FEM and PD domains.

- ✂️ **Distribution logic**: As part of its ongoing development, PeriLab will continue to incorporate new distribution logic for improved performance and reduced computational resources.

- 🏎️ **Optimizations**: As part of its ongoing development, PeriLab will continue to focus on optimizing the simulation process by incorporating new techniques like parallel optimization algorithms for improved efficiency and reduced computational resources.

Feel free to contribute and help us make PeriLab even better! 🙌

## Contributing

We welcome contributions in various forms, including bug reports, documentation improvements, feature suggestions, and more. To get started, follow these steps and have a look at the [Contribution Guidelines](CONTRIBUTING.md):

### Development
1. **Clone the repository:**
```sh
git clone https://gitlab.com/dlr-perihub/PeriLab.jl
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
For any questions or inquiries about PeriLab.jl, feel free to reach out to the authors via email.

## Authors and acknowledgment
[Dr.-Ing. Christian Willberg](mailto::christian.willberg@dlr.de)

[M.Sc. Jan-Timo Hesse](mailto::jan-timo.hesse@dlr.de)

## Project status
`PeriLab` is currently in development.

## Acknowledgments
<p align="center" style="font-size:0;"><!--
  DLR      --><img align="middle" src="https://gitlab.com/dlr-perihub/PeriLab.jl/-/raw/main/assets/dlr.jpg" height="120"><!--
  DFG      --><img align="middle" src="https://gitlab.com/dlr-perihub/PeriLab.jl/-/raw/main/assets/dfg.jpg" height="120"><!--
  SACHSEN  --><img align="middle" src="https://gitlab.com/dlr-perihub/PeriLab.jl/-/raw/main/assets/sachsen.jpg" height="120"><!--
  HyTank  --><img align="middle" src="https://gitlab.com/dlr-perihub/PeriLab.jl/-/raw/main/assets/hytank.jpg" height="120"><!--
  -->
</p>

This project has benefited from funding by the [Deutsche
Forschungsgemeinschaft](https://www.dfg.de/) (DFG, German Research Foundation)
through the following grant ''Gekoppelte Peridynamik-Finite-Elemente-Simulationen zur Schädigungsanalyse von Faserverbundstrukturen''. <br/><br/>Grant number: [WI 4835/5-1](https://gepris.dfg.de/gepris/projekt/456427423)

[M-ERA.NET](https://www.m-era.net/) funded project ''Exploring Multi-Method Analysis of composite structures and joints under consideration of uncertainties engineering and processing (EMMA)''

This measure is co-financed with tax funds on the basis of the budget passed by the [Saxon state parlament](https://www.landtag.sachsen.de/de). <br/><br/>Grant number: [3028223](https://www.m-era.net/materipedia/2020/emma).

[Federal Ministry for Economic Affairs and Climate Action](https://www.bmwk.de/Navigation/DE/Home/home.html) funded project 
''Virtuelle Kennwertermittlung, Schadensprädiktion und Simulationsmethoden für geklebte Fügestellen eines LH2-Tanks in Faserverbundbauweise für die kommerzielle Luftfahrt''.<br/><br/>
Grant number: 20W2214G.