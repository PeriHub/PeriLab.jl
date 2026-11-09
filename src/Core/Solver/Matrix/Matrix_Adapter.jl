# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    Matrix_Adapter

The only PeriLab-aware part of the distributed matrix (`PD_Matrix`): loads the licensed
elements, builds the halo exchange from PeriLab's overlap map, the element data from
`Data_Manager` fields and converts between PeriLab's `nnodes × dof` fields and the
interleaved local vectors of `PD_Matrix` (`x[(i − 1)·dof + c]`).
"""
module Matrix_Adapter

using MPI
using Adapt: adapt
using LinearAlgebra: mul!
using ...Data_Manager
using ...PeriLabExceptions: @abort
using ...ModuleLoader: licensed_modules
using ...ParameterSpec: ParseContext, add_error!, join_path
using ..PD_Matrix: PD_Matrix, HaloExchange, create_halo_exchange, element_type,
                   registered_elements, element_dofs, create_assembly_pattern,
                   create_matrix, create_element_buffer, compute_elements!, assemble!,
                   DistributedOperator, needs_accumulate, select_backend, DofSet,
                   ConstrainedOperator, ShiftedOperator, BiCGStabWorkspace, jacobi!,
                   solve_bicgstab!, constrained_rhs!, apply_prescribed!, exchange!,
                   accumulate_ghosts!, SerialSpace, MPISpace
using ..Correspondence_matrix_based: Correspondence_matrix_based,
                                     precompute_CB_all_nodes!

export load_elements, matrix_element, perilab_halo, correspondence_data, MatrixSystem,
       init_matrix_system, assemble_system!, field_to_vector!, vector_to_field!,
       apply_stiffness!, exchange_field!, update_system!, refresh_active!,
       ImplicitWorkspace, init_implicit, set_prescribed!,
       solve_implicit!

"""
    load_elements()

Includes the licensed matrix elements (license type `MatrixElement`) into
`PD_Matrix.Elements`; they register themselves. Runs once per process.
"""
function load_elements()
    isempty(registered_elements()) || return registered_elements()
    licensed_modules(PD_Matrix.Elements, "MatrixElement")
    return registered_elements()
end

"""
    matrix_element(name::String)

The registered element type `name`; aborts with the list of available elements if the
license does not provide it.
"""
function matrix_element(name::String)
    load_elements()
    T = element_type(name)
    T === nothing &&
        @abort "Matrix element '$name' is not available (licensed feature 'matrix'). Available: $(registered_elements())"
    return T
end

"""
    perilab_halo(dof::Int64; array_type = Array) -> Union{HaloExchange,Nothing}

Halo exchange of the local PeriLab node order (own nodes `1:get_nnodes()`, then the
ghosts), `dof` values per node. The owner of every ghost comes from the overlap map
(`overlap_map[rank + 1][j]["Responder"]`: local ghosts owned by rank `j − 1`). `nothing`
for a single rank.
"""
function perilab_halo(dof::Int64; array_type = Array)
    n_ranks = Data_Manager.get_max_rank()
    n_ranks == 1 && return nothing
    rank = Data_Manager.get_rank()
    nnodes = Data_Manager.get_nnodes()
    n_ghost = Data_Manager.get_num_responder()
    owner = fill(-1, n_ghost)
    for (j, maps) in Data_Manager.get_overlap_map()[rank + 1]
        for local_id in maps["Responder"]
            owner[local_id - nnodes] = j - 1
        end
    end
    any(<(0), owner) && @abort "Ghost nodes without owner in the overlap map."
    return create_halo_exchange(Data_Manager.get_comm(),
                                collect(Int64, Data_Manager.loc_to_glob(1:nnodes)),
                                collect(Int64,
                                        Data_Manager.loc_to_glob((nnodes + 1):(nnodes + n_ghost))),
                                owner, dof; array_type = array_type)
end

"""
    correspondence_data(dof::Int64) -> NamedTuple

Element data of `PeriLab Correspondence` from the `Data_Manager` fields of the own nodes
(flat host arrays, local node numbering; ghosts only appear as neighbours and in
`volume`): the rotated material gradient, the zero-energy stiffness and the `C : B`
tensor of every bond, as in `Correspondence_matrix_based.init_matrix`.
"""
function correspondence_data(dof::Int64)
    nnodes = Data_Manager.get_nnodes()
    nodes = collect(1:nnodes)
    nlist = Data_Manager.get_nlist()
    volume = Data_Manager.get_field("Volume")
    omega = Data_Manager.get_field("Influence Function")
    bond_damage = Data_Manager.get_bond_damage("NP1")
    bond_geometry = Data_Manager.get_field("Bond Geometry")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
    C_voigt = Correspondence_matrix_based._apply_rotation(Data_Manager.get_field("Material Gradient"),
                                                          dof, nnodes)
    use_zero_energy,
    zStiff = Correspondence_matrix_based._setup_zero_energy(nodes, dof, C_voigt,
                                                            inverse_shape_tensor, true)
    CB_nodes = [zeros(length(nlist[i]), dof, dof, dof) for i in nodes]
    precompute_CB_all_nodes!(CB_nodes, nodes, C_voigt, inverse_shape_tensor, nlist,
                             volume, bond_geometry, omega, bond_damage, dof)

    nlist_ptr = cumsum([1; [length(nlist[i]) for i in nodes]])
    n_bonds = nlist_ptr[end] - 1
    geometry = zeros(dof, n_bonds)
    CB = zeros(dof, dof, dof, n_bonds)
    for i in nodes, (k, X) in enumerate(bond_geometry[i])
        b = nlist_ptr[i] + k - 1
        geometry[:, b] .= X
        CB[:, :, :, b] .= CB_nodes[i][k, :, :, :]
    end
    Z = zeros(dof, dof, nnodes)
    use_zero_energy && (Z .= permutedims(zStiff[nodes, :, :], (2, 3, 1)))
    flat(state) = Float64[v for i in nodes for v in state[i]]
    return (nlist_ptr = nlist_ptr, nlist = Int64[v for i in nodes for v in nlist[i]],
            volume = collect(Float64, volume), omega = flat(omega),
            bond_damage = flat(bond_damage), bond_geometry = geometry,
            inverse_shape_tensor = permutedims(inverse_shape_tensor[nodes, :, :],
                                               (2, 3, 1)),
            CB = CB, zStiff = Z), use_zero_energy
end

"""
    MatrixSystem

Distributed stiffness matrix of one rank on its backend: pattern, matrix, element buffer
and data, element, halo exchange, the distributed operator, two local work vectors
(own + ghost degrees of freedom, interleaved) and the status mask (`active` per local
degree of freedom, README A6a.2).
"""
struct MatrixSystem{P,M,V,D,E,H,O,B}
    pattern::P
    K::M
    ke::V
    data::D
    element::E
    halo::H
    op::O
    n_own::Int64
    dof::Int64
    x::V
    y::V
    active::B
end

"""
    init_matrix_system(backend::String) -> MatrixSystem

Symbolic phase on the host, `adapt` to the backend of this rank (`select_backend`), first
assembly. Expects the correspondence fields of the solver setup (material gradient,
inverse shape tensor, bond geometry).
"""
function init_matrix_system(backend::String)
    comm = Data_Manager.get_comm()
    array_type, _, _ = select_backend(backend, comm)
    dof = Data_Manager.get_dof()
    n_own = Data_Manager.get_nnodes()
    n_local = dof * (n_own + Data_Manager.get_num_responder())
    T = matrix_element("PeriLab Correspondence")
    host_data, zero_energy = correspondence_data(dof)
    # invokelatest: the licensed element is defined at runtime (newer world)
    return Base.invokelatest() do
        el = T(dof; zero_energy = zero_energy)
        rows, cols = element_dofs(el, host_data)
        pattern = adapt(array_type, create_assembly_pattern(rows, cols, n_local, n_local))
        halo = perilab_halo(dof; array_type = array_type)
        accumulate = halo !== nothing &&
                     needs_accumulate(comm, pattern.rowptr, dof * n_own)
        K = create_matrix(pattern)
        x = array_type(zeros(n_local))
        sys = MatrixSystem(pattern, K, create_element_buffer(pattern),
                           adapt(array_type, host_data), el, halo,
                           DistributedOperator(K, halo; accumulate = accumulate), n_own,
                           dof, x, similar(x), fill!(similar(x, Bool), true))
        refresh_active!(sys, Data_Manager.get_field("Active"))
        return assemble_system!(sys)
    end
end

"""
    assemble_system!(sys::MatrixSystem)

Numerical phase: element matrices and gather into `sys.K` (pattern unchanged).
"""
function assemble_system!(sys::MatrixSystem)
    Base.invokelatest(compute_elements!, sys.ke, sys.pattern, sys.element, sys.data)
    assemble!(sys.K, sys.pattern, sys.ke; active = sys.active,
              n_own = sys.dof * sys.n_own)
    return sys
end

"""
    field_to_vector!(x, field::AbstractMatrix, n::Int64)

`x[(i − 1)·dof + c] = field[i, c]` for the nodes `1:n` (`dof = size(field, 2)`); device
vectors are filled with one host → device copy.
"""
# the column-major order of the transposed field is the interleaved order
function field_to_vector!(x::Vector{Float64}, field::AbstractMatrix, n::Int64)
    copyto!(x, 1, transpose(view(field, 1:n, :)), 1, n * size(field, 2))
end
function field_to_vector!(x::AbstractVector, field::AbstractMatrix, n::Int64)
    # ponytail: allocates the host copy every call; keep a host buffer if it shows up
    copyto!(x, 1, vec(permutedims(view(field, 1:n, :))), 1, n * size(field, 2))
end

"""
    vector_to_field!(field::AbstractMatrix, x, n::Int64)

Inverse of `field_to_vector!` for the nodes `1:n`.
"""
function vector_to_field!(field::AbstractMatrix, x::Vector{Float64}, n::Int64)
    dof = size(field, 2)
    transpose(view(field, 1:n, :)) .= reshape(view(x, 1:(n * dof)), dof, n)
    return field
end
function vector_to_field!(field::AbstractMatrix, x::AbstractVector, n::Int64)
    return vector_to_field!(field, Array(view(x, 1:(n * size(field, 2)))), n)
end

"""
    apply_stiffness!(f_field, u_field, sys::MatrixSystem)

Internal force of the own nodes `f_field[own, :] = −K u` (exchange of the ghost
displacements, local SpMV, partial sums of the ghost rows to their owners).
"""
function apply_stiffness!(f_field::AbstractMatrix, u_field::AbstractMatrix,
                          sys::MatrixSystem)
    field_to_vector!(sys.x, u_field, sys.n_own)
    mul!(sys.y, sys.op, sys.x)
    vector_to_field!(f_field, sys.y, sys.n_own)
    @views f_field[1:(sys.n_own), :] .*= -1.0
    return f_field
end

"""
    exchange_field!(field::AbstractMatrix, sys::MatrixSystem)

Fills the ghost rows of a `nnodes × dof` field (e.g. the displacements) with the values of
their owners.
"""
function exchange_field!(field::AbstractMatrix, sys::MatrixSystem)
    sys.halo === nothing && return field
    field_to_vector!(sys.x, field, sys.n_own)
    exchange!(sys.halo, sys.x)
    vector_to_field!(field, sys.x, length(sys.x) ÷ sys.dof)
    return field
end

"""
    refresh_active!(sys::MatrixSystem, active::AbstractVector{Bool}) -> Bool

Copies the `Active` field of the own nodes into the status mask, the ghost entries from
their owners. Collective (one `Allreduce`, a halo exchange if anything changed); `true` if
any rank changed, then the matrix has to be reassembled (`update_system!`).
"""
function refresh_active!(sys::MatrixSystem, active::AbstractVector{Bool})
    n, dof = sys.n_own, sys.dof
    host = Array(view(sys.active, 1:(n * dof)))
    changed = any(host[(i - 1) * dof + c] != active[i] for i in 1:n for c in 1:dof)
    sys.halo === nothing || (changed = MPI.Allreduce(changed, |, Data_Manager.get_comm()))
    changed || return false
    values = zeros(length(sys.x))
    for i in 1:n, c in 1:dof
        values[(i - 1) * dof + c] = active[i]
    end
    copyto!(sys.y, values)
    exchange!(sys.halo, sys.y)
    sys.active .= sys.y .> 0.5
    return true
end

"""
    update_system!(sys::MatrixSystem, update::AbstractVector{Bool})

New matrix values after bond damage changed (pattern unchanged, README A6a.1): the element
data (bond damage NP1, inverse shape tensor, `C : B`, zero-energy stiffness) is collected
again from `Data_Manager` and the matrix reassembled. `update` (the `Update` field) tells
whether this rank changed anything, `force` reassembles anyway (status mask changed). No
communication: every rank decides on its own.
"""
function update_system!(sys::MatrixSystem, update::AbstractVector{Bool};
                        force::Bool = false)
    # ponytail: recollects the data of all own nodes; restrict to the flagged nodes (and
    # their neighbours) if damage updates show up in the profile
    force || any(view(update, 1:(sys.n_own))) || return sys
    host, _ = correspondence_data(sys.dof)
    for name in (:bond_damage, :inverse_shape_tensor, :CB, :zStiff)
        copyto!(getfield(sys.data, name), getfield(host, name))
    end
    return assemble_system!(sys)
end

"""
    ImplicitWorkspace

Constrained linear solve with the distributed matrix (`Linear Solver: BiCGStab`):
prescribed mask (local, ghosts from their owners), the operator `K` (or `K + shift`) on the
free degrees of freedom, BiCGStab with Jacobi, `b = −K_fp u_p` of the last solve.
"""
struct ImplicitWorkspace{S<:MatrixSystem,V,B,O,P}
    sys::S
    mask::B
    shift::Union{Nothing,V}
    A::O
    ws::BiCGStabWorkspace{V}
    space::P
    b::V
    rhs::V
    zero::V
end

"""
    init_implicit(sys::MatrixSystem; shift = nothing) -> ImplicitWorkspace

`shift`: per own node (e.g. `a₀ ρ` for Newmark), added to the diagonal of every free
degree of freedom of the node.
"""
function init_implicit(sys::MatrixSystem; shift::Union{Nothing,AbstractVector} = nothing)
    n = sys.dof * sys.n_own
    mask = fill!(similar(sys.x, Bool), false)
    prescribed = DofSet(mask)
    shift_vector = nothing
    if shift !== nothing
        host = zeros(length(sys.x))
        for i in 1:(sys.n_own), c in 1:(sys.dof)
            host[(i - 1) * sys.dof + c] = shift[i]
        end
        shift_vector = copyto!(similar(sys.x), host)
    end
    op = ConstrainedOperator(sys.K, prescribed; n_own = n)
    op = shift_vector === nothing ? op : ShiftedOperator(op, shift_vector)
    A = DistributedOperator(op, sys.halo; accumulate = sys.op.accumulate)
    space = sys.halo === nothing ? SerialSpace(n) : MPISpace(Data_Manager.get_comm(), n)
    return ImplicitWorkspace(sys, mask, shift_vector, A, BiCGStabWorkspace(sys.x), space,
                             similar(sys.x), similar(sys.x), fill!(similar(sys.x), 0.0))
end

"""
    set_prescribed!(iw::ImplicitWorkspace, free::AbstractVector{Int64}, active::AbstractVector{Bool})

Prescribed degrees of freedom: all except `free` (`Data_Manager.get_bc_free_dof()`,
component-major `i + (c − 1)·nnodes`), plus every degree of freedom of inactive own nodes;
the ghost entries come from their owners. Recomputes the Jacobi preconditioner.
"""
function set_prescribed!(iw::ImplicitWorkspace, free::AbstractVector{Int64},
                         active::AbstractVector{Bool})
    sys = iw.sys
    n, dof = sys.n_own, sys.dof
    host = ones(length(sys.x))
    for g in free
        i, c = mod1(g, n), (g - 1) ÷ n + 1
        host[(i - 1) * dof + c] = 0.0
    end
    for i in 1:n, c in 1:dof
        active[i] || (host[(i - 1) * dof + c] = 1.0)
    end
    copyto!(sys.y, host)
    exchange!(sys.halo, sys.y)
    iw.mask .= sys.y .> 0.5
    jacobi!(iw.ws, sys.K, DofSet(iw.mask); halo = sys.halo,
            accumulate = sys.op.accumulate, n_own = n * dof, shift = iw.shift)
    return iw
end

"""
    solve_implicit!(u, f, iw::ImplicitWorkspace; rtol, itmax, b_field = nothing)

Solves `(K + shift) u_f = f_f − K_fp u_p` for the free degrees of freedom of the own nodes
with BiCGStab, starting from the current `u_f`; `u` holds the prescribed values. `b_field`
receives `−K_fp u_p` (free entries, prescribed entries zero).

# Returns
- `stats::NamedTuple`: `solved`, `niter`, `residual` (relative)
"""
function solve_implicit!(u::AbstractMatrix, f::AbstractMatrix, iw::ImplicitWorkspace;
                         rtol::Float64, itmax::Int64,
                         b_field::Union{Nothing,AbstractMatrix} = nothing)
    sys = iw.sys
    n = sys.n_own
    prescribed = DofSet(iw.mask)
    field_to_vector!(sys.x, u, n)
    exchange!(sys.halo, sys.x)
    constrained_rhs!(iw.b, sys.K, prescribed, sys.x, iw.zero)
    sys.op.accumulate && accumulate_ghosts!(sys.halo, iw.b)
    field_to_vector!(iw.rhs, f, n)
    @. iw.rhs = ifelse(iw.mask, 0.0, iw.rhs + iw.b)
    # start vector: current free values, zero at the prescribed degrees of freedom
    @. sys.y = ifelse(iw.mask, 0.0, sys.x)
    x, stats = solve_bicgstab!(iw.ws, iw.A, iw.rhs, iw.space; rtol = rtol, itmax = itmax,
                               x0 = sys.y)
    apply_prescribed!(x, prescribed, sys.x)
    vector_to_field!(u, x, n)
    if b_field !== nothing
        @. iw.b = ifelse(iw.mask, 0.0, iw.b)
        vector_to_field!(b_field, iw.b, n)
    end
    return stats
end

"Input checks shared by the implicit matrix solvers (Linear Static Matrix Based, Newmark)."
function check_implicit_backend!(p, path::String, ctx::ParseContext)
    distributed = p.matrix_backend != "SparseArrays"
    distributed == (p.linear_solver == "BiCGStab") ||
        add_error!(ctx, join_path(path, "Linear Solver"),
                   "Linear Solver: BiCGStab requires Matrix Backend CPU, GPU or Mixed, these backends require Linear Solver: BiCGStab")
    distributed && p.matrix_update &&
        add_error!(ctx, join_path(path, "Matrix Update"),
                   "Matrix Update requires Matrix Backend: SparseArrays")
    return nothing
end

end
