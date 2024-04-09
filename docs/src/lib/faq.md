# PeriLab - Frequently Asked Questions (FAQ)

## General Questions

### What is PeriLab?

PeriLab is a peridynamic simulation software designed for different kind of mechanical problems.

## Installation and Setup

### What are the system requirements for installing PeriLab?

In order to install PeriLab, you will need to have a recent version of Julia or Docker installed on your system.

### How do I update my PeriLab software?

To update PeriLab you just need to run `julia -e "using Pkg; Pkg.update()"` or pull the latest version of the PeriLab Docker image.

## Simulation and Analysis

### How can I generate my own input mesh?

There are multiple options to generate your own input mesh:

- Use the [Perihub](https://github.com/PeriHub/Perihub) framework to generate your own input mesh.

- Generate your own input mesh with julia, [here](https://github.com/PeriHub/PeriLab.jl/blob/main/examples/Training/meshing/mesh.jl) is an example.

- Create a .png image of your 2D model and translate it with this [script](https://github.com/PeriHub/PeriLab.jl/blob/main/examples/Training/meshing/image.jl).

- Use an existing Abaqus Inputfile (.inp).

- Transfer your mesh using [meshio](https://github.com/nschloe/meshio)

- Create a feature issue and let us know what you need.

### What to do with my results?

First of all congratulations 🎉 on the success of your simulation. Now, you can take a look at your results. To do that, you can use [ParaView](https://www.paraview.org/) it's free and very powerful.

### Can I customize parameters for my simulations in PeriLab?

Yes, PeriLab allows users to customize various parameters to tailor simulations based on their specific requirements.

## Troubleshooting

### I'm experiencing technical issues with PeriLab. What should I do?

If you encounter technical issues, please create an issue and describe it in detail.

## Additional Assistance

If your question is not addressed here, please feel free to [contact us](mailto:christian.willberg@dlr.de) or create an issue for further assistance.

## Contributions

Contributions are always welcomed, take a look at the [Contributing Guidelines](https://github.com/PeriHub/PeriLab.jl/blob/main/CONTRIBUTING.md)
