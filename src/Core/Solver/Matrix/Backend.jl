# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Backend

using MPI

export select_backend, GPU_HOOK, rank_weights

"""
    GPU_HOOK

`nothing` until the CUDA extension (`ext/PeriLabCUDAExt.jl`) is loaded; then a
`NamedTuple` `(array_type, n_devices, device!)`.
"""
const GPU_HOOK = Ref{Any}(nothing)

const CUDA_ID = Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA")

"Loads CUDA.jl (and with it the extension) if it is installed but not loaded yet."
function load_cuda!()
    GPU_HOOK[] === nothing || return nothing
    try
        Base.require(CUDA_ID)
    catch
        # not installed: GPU_HOOK stays nothing
    end
    return nothing
end

"""
    select_backend(mode::String, comm::MPI.Comm) -> (array_type, on_gpu, device)

Array type of this rank for `mode` (`"CPU"`, `"GPU"`, `"Mixed"`) and its CUDA device
(`-1` on the CPU). GPUs are assigned by the rank within the node:
`device = local_rank % number_of_devices`. `"GPU"` without a functional GPU throws;
`"Mixed"` puts one rank per GPU on the GPU and the other ranks on the CPU, without a GPU
every rank runs on the CPU (with one warning).
"""
function select_backend(mode::String, comm::MPI.Comm)
    mode in ("CPU", "GPU", "Mixed") ||
        throw(ArgumentError("Unknown Matrix Backend $mode (CPU | GPU | Mixed)."))
    mode == "CPU" && return (Array, false, -1)
    load_cuda!()
    hook = GPU_HOOK[]
    # invokelatest: the hook functions are defined by the extension at runtime
    n_devices = hook === nothing ? 0 : Base.invokelatest(hook.n_devices)
    if n_devices == 0
        mode == "GPU" &&
            error("Matrix Backend GPU requires CUDA.jl (using CUDA) and a functional GPU")
        MPI.Comm_rank(comm) == 0 &&
            @warn "Matrix Backend Mixed without a functional GPU (using CUDA): all ranks run on the CPU"
        return (Array, false, -1)
    end
    rank = MPI.Comm_rank(comm)
    node_comm = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, rank)
    local_rank = MPI.Comm_rank(node_comm)
    if mode == "GPU" || local_rank < n_devices
        device = local_rank % n_devices
        Base.invokelatest(hook.device!, device)
        return (hook.array_type, true, device)
    end
    return (Array, false, -1)
end

"""
    rank_weights(comm::MPI.Comm, backend::String, gpu_weight::Float64) -> Vector{Float64}

Relative number of nodes of every rank: `gpu_weight` for the ranks `select_backend` puts
on a GPU, 1 otherwise (one `Allgather`).
"""
function rank_weights(comm::MPI.Comm, backend::String, gpu_weight::Float64)
    _, on_gpu, _ = select_backend(backend, comm)
    return MPI.Allgather(on_gpu ? gpu_weight : 1.0, comm)
end

end
