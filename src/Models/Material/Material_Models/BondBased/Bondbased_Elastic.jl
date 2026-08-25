# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Elastic

using .......Data_Manager
using .......PeriLabExceptions: @abort
using ....Material_Basis: get_symmetry, apply_pointwise_E, compute_bond_based_constants
using .......Helpers: is_dependent
using LoopVectorization
using TimerOutputs: @timeit
using KernelAbstractions
using Adapt
using CUDA
export init_model
export fe_support
export material_name
export compute_model

"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
    return false
end

"""
  init_model(nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
"""
function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    constant = Data_Manager.create_constant_node_scalar_field("Bond Based Constant",
                                                              Float64)
    horizon = Data_Manager.get_field("Horizon")
    symmetry::String = get_symmetry(material_parameter)
    compute_bond_based_constants(nodes, symmetry, constant, horizon)
end

"""
	material_name()

Returns the name of the material model.
"""
function material_name()
    return "Bond-based Elastic"
end

"""
    select_backend(name::AbstractString)

Resolves the `"Backend"` material parameter to a concrete
`KernelAbstractions` backend.

- `"CPU"` (default if the parameter is not set) -> `KernelAbstractions.CPU()`
- `"CUDA"` -> `CUDA.CUDABackend()`, provided `CUDA.functional()` reports a
  usable GPU on this machine.

`CUDA.jl` is imported unconditionally at the top of this module. This is
safe even on machines without an NVIDIA GPU or driver -- `using CUDA` does
not error in that case, it simply leaves `CUDA.functional()` returning
`false`, which is what this function checks before ever touching
`CUDABackend()`.

Any unrecognized or currently-unusable choice aborts with a clear message
via PeriLab's own `@abort` convention, rather than failing later with an
opaque `MethodError` from inside a kernel launch.
"""
function select_backend(name::AbstractString)
    name_lc = lowercase(name)
    if name_lc == "cpu"
        return KernelAbstractions.CPU()
    elseif name_lc == "cuda"
        if !CUDA.functional()
            @abort "Backend \"CUDA\" was requested for Bondbased_Elastic, but CUDA.jl reports no functional GPU on this system (check drivers/toolkit installation). Set \"Backend\" to \"CPU\" or fix the GPU setup."
            return KernelAbstractions.CPU()
        end
        return CUDA.CUDABackend()
    else
        @abort "Unknown \"Backend\" option \"$name\" for Bondbased_Elastic. Use \"CPU\" or \"CUDA\"."
        return KernelAbstractions.CPU()
    end
end

"""
	compute_model(nodes::AbstractVector{Int64}, material_parameter::Dict, time::Float64, dt::Float64)

Calculate the elastic bond force for each node.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.

The `material_parameter` dict may include an optional `"Backend"` entry
(`"CPU"` or `"CUDA"`, case-insensitive) to select where the bond force
kernel runs. Defaults to `"CPU"` if not given, which reproduces the
original CPU-only behavior (aside from the flatten/scatter overhead noted
in `compute_bb_force!`'s docstring).
"""
function compute_model(nodes::AbstractVector{Int64},
                       material_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    constant::NodeScalarField{Float64} = Data_Manager.get_field("Bond Based Constant")

    undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
    deformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Deformed Bond Geometry",
                                                                     "NP1")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")
    bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    bond_force::BondVectorState{Float64} = Data_Manager.get_field("Bond Forces")

    E = material_parameter["Young's Modulus"]

    dependend_value,
    dependent_field = is_dependent("Young's Modulus", material_parameter)

    @timeit "any zero" begin
        for iID in nodes
            for deformed_length in deformed_bond_length[iID]
                if deformed_length == 0
                    @abort "Length of bond is zero due to its deformation."
                    return nothing
                end
            end
        end
    end

    backend = select_backend(get(material_parameter, "Backend", "CPU"))

    # Calculate the bond force for every node in `nodes` in a single kernel
    # launch (see compute_bb_force! docstring for how the ragged per-node
    # data is flattened for this).
    @timeit "compute_bb_force!" compute_bb_force!(nodes, bond_force, constant,
                                                  bond_damage, deformed_bond_length,
                                                  undeformed_bond_length, deformed_bond;
                                                  backend = backend)

    # might be put in constant
    if dependend_value
        @timeit "apply_pointwise_E" apply_pointwise_E(nodes, E, bond_force, dependent_field)
    else
        @timeit "apply_pointwise_E" apply_pointwise_E(nodes, E, bond_force)
    end
end

"""
    compute_bb_force!(nodes, bond_force, constant, bond_damage, deformed_bond_length,
                       undeformed_bond_length, deformed_bond; backend=KernelAbstractions.CPU())

Calculates the bond-based elastic force for every node in `nodes`, using a
single `KernelAbstractions` kernel so the SAME implementation runs
unchanged on CPU or CUDA -- pass `backend=KernelAbstractions.CPU()` (the
default) or `backend=CUDA.CUDABackend()`.

# Why flatten
`bond_force`, `bond_damage`, `deformed_bond_length`, `undeformed_bond_length`,
and `deformed_bond` are all stored per-node as ragged/jagged arrays (each
node has a different number of neighbors). GPU kernels need flat,
fixed-stride arrays, so this function:

  1. Flattens the relevant per-node data for `nodes` into plain arrays --
     one entry per (node, neighbor) bond -- keeping parallel `owner_node`/
     `owner_neighbor` index arrays so results can be scattered back later.
  2. Moves those flat arrays onto the target `backend` via `Adapt.adapt`
     (a no-op on `CPU()`, an actual host->device copy on `CUDABackend()`).
  3. Launches `bb_force_kernel!` ONCE over all bonds at once -- one
     GPU/CPU thread per bond. Each thread only ever writes to its own
     unique output column, so this is race-free with no atomics needed.
  4. Copies the result back to the host (a no-op on CPU) and scatters it
     into the original ragged `bond_force` structure so the rest of
     PeriLab (which still expects that layout) is unaffected.

# Performance note
Steps 1 and 4 run on the CPU every call and are new overhead that the
original per-node-loop implementation did not have -- this is the cost of
keeping the change contained to this one file/function rather than
restructuring PeriLab's Data_Manager to store bonds in a flat/CSR layout
permanently (the natural follow-up, which would let this flatten happen
once when the neighbor list is built instead of every force evaluation).
On `backend=CPU()` this is a straightforward trade of some setup cost for
having one unified, GPU-capable implementation instead of two separate
code paths to maintain. On `backend=CUDA.CUDABackend()`, this overhead is
expected to be increasingly outweighed by GPU throughput as problem size
grows; it has NOT been benchmarked as part of this change (per request,
this edit was not run/tested).
"""
function compute_bb_force!(nodes::AbstractVector{Int64},
                           bond_force::BondVectorState{Float64},
                           constant::NodeScalarField{Float64},
                           bond_damage::BondScalarState{Float64},
                           deformed_bond_length::BondScalarState{Float64},
                           undeformed_bond_length::BondScalarState{Float64},
                           deformed_bond::BondVectorState{Float64};
                           backend = KernelAbstractions.CPU())
    total_bonds = 0
    dof = 0
    @inbounds for iID in nodes
        n_i = length(bond_damage[iID])
        total_bonds += n_i
        if dof == 0 && n_i > 0
            dof = length(deformed_bond[iID][1])
        end
    end
    total_bonds == 0 && return nothing

    flat_constant = Vector{Float64}(undef, total_bonds)
    flat_damage = Vector{Float64}(undef, total_bonds)
    flat_deformed_length = Vector{Float64}(undef, total_bonds)
    flat_undeformed_length = Vector{Float64}(undef, total_bonds)
    flat_deformed_bond = Matrix{Float64}(undef, dof, total_bonds)
    # bookkeeping needed to scatter kernel results back into the ragged
    # bond_force structure afterwards
    owner_node = Vector{Int64}(undef, total_bonds)
    owner_neighbor = Vector{Int64}(undef, total_bonds)

    idx = 0
    @inbounds for iID in nodes
        n_i = length(bond_damage[iID])
        for k in 1:n_i
            idx += 1
            flat_constant[idx] = 0.5 * constant[iID]
            flat_damage[idx] = bond_damage[iID][k]
            flat_deformed_length[idx] = deformed_bond_length[iID][k]
            flat_undeformed_length[idx] = undeformed_bond_length[iID][k]
            for d in 1:dof
                flat_deformed_bond[d, idx] = deformed_bond[iID][k][d]
            end
            owner_node[idx] = iID
            owner_neighbor[idx] = k
        end
    end

    # Move to the target device. Identity (no copy) when backend is CPU().
    dev_constant = Adapt.adapt(backend, flat_constant)
    dev_damage = Adapt.adapt(backend, flat_damage)
    dev_deformed_length = Adapt.adapt(backend, flat_deformed_length)
    dev_undeformed_length = Adapt.adapt(backend, flat_undeformed_length)
    dev_deformed_bond = Adapt.adapt(backend, flat_deformed_bond)
    dev_bond_force = Adapt.adapt(backend, similar(flat_deformed_bond))

    kernel! = bb_force_kernel!(backend)
    kernel!(dev_bond_force, dev_constant, dev_damage, dev_deformed_length,
            dev_undeformed_length, dev_deformed_bond, dof; ndrange = total_bonds)
    KernelAbstractions.synchronize(backend)

    # Bring the result back to the host (no-op on CPU) and scatter it into
    # the ragged bond_force structure the rest of PeriLab expects.
    host_bond_force = Array(dev_bond_force)
    @inbounds for idx in 1:total_bonds
        iID = owner_node[idx]
        k = owner_neighbor[idx]
        for d in 1:dof
            bond_force[iID][k][d] = host_bond_force[d, idx]
        end
    end

    return nothing
end

"""
    bb_force_kernel!(bond_force, constant, damage, deformed_length,
                      undeformed_length, deformed_bond, dof)

`KernelAbstractions` kernel computing the classic bond-based elastic (PMB)
force for one bond per thread:

    scalar = constant * damage * (deformed_length - undeformed_length) / undeformed_length
    bond_force[:, idx] = scalar * deformed_bond[:, idx] / deformed_length

`constant` here already includes the caller's `0.5 *` factor, matching the
original (pre-GPU) `compute_bb_force!` implementation. Every thread writes
only to its own column of `bond_force`, so no atomics are required.
"""
@kernel function bb_force_kernel!(bond_force, @Const(constant), @Const(damage),
                                  @Const(deformed_length), @Const(undeformed_length),
                                  @Const(deformed_bond), dof::Int)
    idx = @index(Global)

    scalar = constant[idx] * damage[idx] *
             (deformed_length[idx] - undeformed_length[idx]) / undeformed_length[idx]

    @inbounds for d in 1:dof
        bond_force[d, idx] = scalar * deformed_bond[d, idx] / deformed_length[idx]
    end
end

"""
	fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.

# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
