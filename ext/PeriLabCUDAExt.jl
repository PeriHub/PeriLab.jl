# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    PeriLabCUDAExt

GPU backend of the distributed stiffness matrix (`Matrix Backend: GPU / Mixed`): registers
`CuArray` and the device selection in `PD_Matrix.GPU_HOOK`.
"""
module PeriLabCUDAExt

using CUDA: CUDA, CuArray
using PeriLab

function __init__()
    PeriLab.Solver_Manager.PD_Matrix.GPU_HOOK[] = (array_type = CuArray,
                                                   n_devices = () -> CUDA.functional() ?
                                                                     length(CUDA.devices()) :
                                                                     0,
                                                   device! = d -> (CUDA.device!(d); nothing
                                                                   ))
    return nothing
end

end
