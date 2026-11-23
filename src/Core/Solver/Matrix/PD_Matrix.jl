# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    PD_Matrix

Sparse, non-symmetric peridynamic system matrices on CPU and GPU (KernelAbstractions),
distributed over MPI with one-horizon ghosts. No access to PeriLab modules: all data comes
through function arguments (prototype README A0.3); `Matrix_Adapter` connects it to
`Data_Manager`. The element implementations are licensed modules that register
themselves with `register_element` (see `Models/Elements/Element_template`).
"""
module PD_Matrix

include("Core/Sparse_CSR.jl")
include("Core/Matrix_Assembly.jl")
include("Core/Dof_Sets.jl")
include("Core/Operators.jl")
include("MPI_communication/Halo_Exchange.jl")
include("Core/Distributed_Operators.jl")
include("Core/Solver/BiCGStab.jl")
include("Models/Elements/Elements.jl")
include("Backend.jl")

using .Sparse_CSR
using .Matrix_Assembly
using .Dof_Sets
using .Operators
using .Halo_Exchange
using .Distributed_Operators
using .BiCGStab
using .Elements
using .Backend

export CSRMatrix, to_sparse
export AssemblyPattern, create_assembly_pattern, create_matrix, create_element_buffer,
       assemble!
export DofSet, create_dof_set
export ConstrainedOperator, ShiftedOperator, constrained_rhs!, apply_prescribed!
export HaloExchange, create_halo_exchange, exchange!, accumulate_ghosts!
export AbstractSpace, SerialSpace, MPISpace, inner, DistributedOperator, needs_accumulate
export BiCGStabWorkspace, jacobi!, solve_bicgstab!
export AbstractElement, element_name, element_dofs, compute_element!, compute_elements!,
       n_elements, register_element, element_type, registered_elements
export select_backend, GPU_HOOK, rank_weights

"""
    build_system(el, data, n_rows, n_cols = n_rows)

Convenience for the symbolic phase on the host: pattern, matrix and element buffer for
element `el` with host data `data`.
"""
function build_system(el::AbstractElement, data, n_rows::Int64, n_cols::Int64 = n_rows)
    rows, cols = element_dofs(el, data)
    pattern = create_assembly_pattern(rows, cols, n_rows, n_cols)
    return pattern, create_matrix(pattern), create_element_buffer(pattern)
end
export build_system

end
