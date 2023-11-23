# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2

module use /home/f_peridi/julia/easybuild-installed/modules/all
module load Julia/1.9.3-linux-x86_64
module load OpenMPI/3.1.3

srun julia --project=. src/main.jl examples/Undulations/Smetana.yaml -v