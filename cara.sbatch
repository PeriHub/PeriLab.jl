#!/bin/bash

# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#SBATCH --time=01:00:00
#SBATCH --ntasks=50

module use /home/f_peridi/julia/easybuild-installed/modules/all
module load Julia/1.9.3-linux-x86_64
module load OpenMPI/3.1.3

srun julia --project=. src/main.jl examples/Benchmark/Benchmark_3dx_1000.yaml -v
srun julia --project=. src/main.jl examples/Benchmark/Benchmark_4dx_1000.yaml -v
srun julia --project=. src/main.jl examples/Benchmark/Benchmark_5dx_1000.yaml -v