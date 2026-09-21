# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Serialization

export write_checkpoint
export read_checkpoint!

# Keys in `data` that must NOT be restored from disk, because they hold
# things that are only valid for the process that created them (an MPI
# communicator, live Module references) or are re-derived every run.
const CHECKPOINT_EXCLUDED_KEYS = ("commMPi", "model_modules", "active_models",
                                  "all_active_models", "rank", "max_rank",
                                  "mpi_active")

"""
    checkpoint_file(directory::String, rank::Int64)

Path of the restart file for a given MPI rank. One file per rank, so any
future `mpiexecjl` launch with the same rank count can restore its own
partition independently.
"""
function checkpoint_file(directory::String, rank::Int64)
    return joinpath(directory, "restart", "checkpoint_rank$(rank).jls")
end

"""
    write_checkpoint(directory::String, rank::Int64)

Serializes this rank's field arrays (`fieldmanager.fields`) and the
persistable part of `data` (current time, step, field bookkeeping, etc.) to
disk. Intended to be called between solver steps, so a fresh
`mpiexecjl` invocation started with `--reload` can pick up where this one
left off.
"""
function write_checkpoint(directory::String, rank::Int64)
    path = checkpoint_file(directory, rank)
    mkpath(dirname(path))
    state = Dict(k => v for (k, v) in data if !(k in CHECKPOINT_EXCLUDED_KEYS))
    tmp_path = path * ".tmp"
    open(tmp_path, "w") do io
        serialize(io, (fields = fieldmanager.fields, state = state))
    end
    mv(tmp_path, path; force = true)
end

"""
    read_checkpoint!(directory::String, rank::Int64)

Restores this rank's `data` and `fieldmanager.fields` from a checkpoint
written by `write_checkpoint`. Always runs `initialize_data()` first so
excluded/derived keys (rank, comm, model wiring) start clean and get
re-populated normally later in `run()`.

Falls back to a plain `initialize_data()` (i.e. a fresh start) if no
checkpoint file is found, so passing `--reload` on a first run is harmless.
"""
function read_checkpoint!(directory::String, rank::Int64)
    path = checkpoint_file(directory, rank)
    if !isfile(path)
        @info "No checkpoint found; starting fresh instead of reloading."
        return
    end
    checkpoint = open(deserialize, path)
    for (k, v) in checkpoint.state
        data[k] = v
    end
    empty!(fieldmanager.fields)
    merge!(fieldmanager.fields, checkpoint.fields)
end
