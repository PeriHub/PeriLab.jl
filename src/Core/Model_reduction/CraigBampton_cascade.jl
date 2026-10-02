# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module CraigBampton_cascade
using LinearAlgebra
using SparseArrays
using TimerOutputs: @timeit
export model_reduction_name
export reduce_matrices

"""
    model_reduction_name()

Name of the reduction scheme, used to identify it in the input deck and in log output.

# Returns
- `::String`: The name of the scheme
"""
function model_reduction_name()
    return "Craig Bampton Cascade"
end

"""
    nonzero_rows(A)

Indices of the rows of a matrix that hold at least one entry.

A condensed degree of freedom whose row is empty carries no stiffness at all — an
isolated node, or one whose bonds were all cut by a bond filter. Such a row makes `Kcc`
singular, and the eigenvalue problem cannot separate its rigid body mode from the modes
being looked for. They are removed before the factorization and the eigenvalue problem,
and the corresponding entries of the modes stay zero.

# Arguments
- `A::AbstractMatrix`: The matrix
# Returns
- `::Vector{Int64}`: Indices of the non-empty rows
"""
function nonzero_rows(A::SparseMatrixCSC)
    occupied = falses(size(A, 1))
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for column in axes(A, 2)
        for index in nzrange(A, column)
            if values[index] != 0.0
                occupied[rows[index]] = true
            end
        end
    end
    return findall(occupied)
end

nonzero_rows(A::AbstractMatrix) = findall(!iszero, vec(sum(abs, A; dims = 2)))

"""
    nonzero_columns(A)

Indices of the columns of a matrix that hold at least one entry.

For `Kcr` these are the retained degrees of freedom with a bond into the condensed
region — the coupling layer. Every other column is empty, so its recovery mode is zero
and never has to be solved for. That is the difference between one triangular solve per
retained degree of freedom and one per coupling degree of freedom, and it also bounds
where the reduced matrices can fill in: `Krc * B` has entries only in the rows and
columns of the coupling layer, everything else keeps the sparsity of `Krr`.

# Arguments
- `A::AbstractMatrix`: The matrix
# Returns
- `::Vector{Int64}`: Indices of the non-empty columns
"""
function nonzero_columns(A::SparseMatrixCSC)
    values = nonzeros(A)
    columns = Int64[]
    for column in axes(A, 2)
        for index in nzrange(A, column)
            if values[index] != 0.0
                push!(columns, column)
                break
            end
        end
    end
    return columns
end

nonzero_columns(A::AbstractMatrix) = findall(!iszero, vec(sum(abs, A; dims = 1)))

"""
    _warn_not_pd(err = nothing)

Warns that the condensed stiffness failed to factorize as positive definite, shared
between the dense in-place path and the sparse fallback below.
"""
function _warn_not_pd(err = nothing)
    @warn "Cholesky of the condensed stiffness failed, falling back to LU. Kcc is not " *
          "positive definite: check whether every condensed node is tied to the " *
          "retained region -- an isolated node or a detached cluster has a rigid body " *
          "mode." exception=err
end

"""
    factorize_condensed!(Kcc)

Factorization of the condensed stiffness, Cholesky where possible.

With the interface held fixed `Kcc` is symmetric positive definite. A sparse Cholesky is
faster to build than an LU and considerably faster to apply to a dense block of right
hand sides, because CHOLMOD goes through supernodal BLAS3 kernels while UMFPACK works
column by column. With one right hand side per coupling degree of freedom that
difference dominates the whole reduction.

A failure is informative in itself: `Kcc` is then not positive definite, which for a
fixed interface points at a condensed region falling apart into pieces not tied to any
retained degree of freedom — or, in peridynamics, simply at a stiffness that is not
symmetric.

Dense `Kcc` (a `StridedMatrix` of a BLAS float type -- what a single, small condensed
region typically comes down to) factorizes in place: LAPACK's `potrf` touches only the
upper triangle it is asked for, so on failure the lower triangle plus the saved diagonal
is enough to rebuild the original matrix for the LU fallback, without a second copy of
`Kcc`. Anything else -- in particular a `SparseMatrixCSC`, which `cholesky!`/`lu!` do not
support in place -- goes through the ordinary, allocating `cholesky`/`lu`.

# Arguments
- `Kcc::AbstractMatrix`: Stiffness of the condensed part
# Returns
- A factorization object supporting `ldiv!`
"""
function factorize_condensed!(Kcc::StridedMatrix{<:LinearAlgebra.BlasFloat})
    n = LinearAlgebra.checksquare(Kcc)
    d = diag(Kcc)                                   # nur die Diagonale sichern
    F = cholesky!(Symmetric(Kcc, :U); check = false)
    issuccess(F) && return F

    _warn_not_pd()
    # potrf hat nur das obere Dreieck inkl. Diagonale verändert; das untere ist intakt
    @inbounds for j in 1:n
        Kcc[j, j] = d[j]
        for i in 1:(j - 1)
            Kcc[i, j] = Kcc[j, i]
        end
    end
    return lu!(Kcc)
end

function factorize_condensed!(Kcc::AbstractMatrix)
    try
        return cholesky(Symmetric(Kcc))
    catch err
        _warn_not_pd(err)
        return lu(Kcc)
    end
end

"""
    check_symmetry(K; samples = 2000, tolerance = 1e-8)

Measures how far the stiffness is from symmetric, sampling stored entries.

Craig-Bampton is built on symmetry: the fixed-interface modes use the symmetric solver,
and the classical derivation drops the stiffness coupling blocks because `Kcc B = Kcr`
makes them cancel — which needs `Krc = Kcr'`. In peridynamics the stiffness is not
symmetric in general, two points with different horizons do not contribute equally to
each other, so the deviation is measured rather than assumed away. The coupling blocks
are formed explicitly either way; this only reports how much they will carry.

Sampling runs over the stored entries of `K`, not over random index pairs. For a sparse
matrix a random pair is almost surely a structural zero on both sides, which would make
any deviation invisible and report perfect symmetry for every input.

# Arguments
- `K::AbstractMatrix`: The full stiffness matrix
# Keywords
- `samples::Int64`: Number of stored entries drawn
- `tolerance::Float64`: Relative deviation accepted as rounding
# Returns
- `::Float64`: The largest relative deviation found
"""
function check_symmetry(K::SparseMatrixCSC; samples::Int64 = 2000,
                        tolerance::Float64 = 1e-8)
    rows = rowvals(K)
    values = nonzeros(K)
    total = length(values)
    total == 0 && return 0.0

    stride = max(1, total ÷ samples)
    scale = 0.0
    deviation = 0.0
    drawn = 0

    @inbounds for column in axes(K, 2)
        for index in nzrange(K, column)
            index % stride == 0 || continue
            i = rows[index]
            upper = values[index]
            lower = K[column, i]
            scale = max(scale, abs(upper), abs(lower))
            deviation = max(deviation, abs(upper - lower))
            drawn += 1
        end
    end

    drawn == 0 && return 0.0
    relative = deviation / max(scale, eps())

    if relative > tolerance
        @warn "Stiffness matrix is not symmetric (largest relative deviation " *
              "$(round(relative; sigdigits = 3)) over $drawn sampled entries). The " *
              "modal coupling blocks are kept, so the reduction stays consistent, but " *
              "the fixed-interface modes are those of the symmetric part."
    end
    return relative
end

check_symmetry(K::AbstractMatrix; kwargs...) = check_symmetry(sparse(K); kwargs...)

"""
    report_modal_coupling(K_reduced, M_reduced, nr, n_modes)

Reports the coupling between the retained and the modal part, and warns if it is empty.

A mode that no block couples to stays at rest for the whole simulation. It still carries
mass and shifts the eigenfrequencies, but never responds, which makes the result worse
than the static condensation it was meant to improve on — and, characteristically,
independent of the number of modes.

There are two paths. `Mrcm` drives the modes through the acceleration of the retained
degrees of freedom and is always present. `Krcm` drives them through the displacement and
vanishes exactly for a symmetric stiffness, so an empty one is only expected there.

# Arguments
- `K_reduced::SparseMatrixCSC`, `M_reduced::SparseMatrixCSC`: The reduced matrices
- `nr::Int64`: Number of retained degrees of freedom
- `n_modes::Int64`: Number of modes
# Returns
- `::Bool`: Whether any coupling holds entries
"""
function report_modal_coupling(K_reduced::SparseMatrixCSC, M_reduced::SparseMatrixCSC,
                               nr::Int64, n_modes::Int64)
    n_modes == 0 && return true

    modal = (nr + 1):(nr + n_modes)
    mass_coupling = M_reduced[1:nr, modal]
    stiffness_coupling = K_reduced[1:nr, modal]

    mass_entries = nnz(mass_coupling)
    stiffness_entries = nnz(stiffness_coupling)

    if mass_entries == 0 && stiffness_entries == 0
        @warn "Neither the mass nor the stiffness couples the retained and the modal " *
              "part. The modes will stay at rest and the reduction behaves worse than " *
              "a static condensation."
        return false
    end

    retained_mass = maximum(abs, nonzeros(M_reduced[1:nr, 1:nr]); init = 0.0)
    retained_stiffness = maximum(abs, nonzeros(K_reduced[1:nr, 1:nr]); init = 0.0)

    @info "Modal coupling: mass $mass_entries entries, largest " *
          "$(round(maximum(abs, nonzeros(mass_coupling); init = 0.0); sigdigits = 3)) " *
          "against $(round(retained_mass; sigdigits = 3)); stiffness " *
          "$stiffness_entries entries, largest " *
          "$(round(maximum(abs, nonzeros(stiffness_coupling); init = 0.0); sigdigits = 3)) " *
          "against $(round(retained_stiffness; sigdigits = 3))"

    return true
end

"""
    graph_distance(K0, nodes, dof, local_index, sources)

Breadth first search over the graph of `K0` restricted to `nodes`, starting from the
local indices `sources` at distance 1.

# Returns
- `::Vector{Int64}`: Distance of every node, `typemax(Int64)` where unreachable
"""
function graph_distance(K0::SparseMatrixCSC{Float64,Int64}, nodes::Vector{Int64},
                        dof::Int64, local_index::Vector{Int64}, sources::Vector{Int64})
    n_all = size(K0, 1) ÷ dof
    rows = rowvals(K0)
    values = nonzeros(K0)
    distance = fill(typemax(Int64), length(nodes))
    distance[sources] .= 1
    frontier = copy(sources)
    level = 1
    while !isempty(frontier)
        level += 1
        next = Int64[]
        for j in frontier, d in 1:dof
            column = (d - 1) * n_all + nodes[j]
            for index in nzrange(K0, column)
                values[index] == 0.0 && continue
                k = local_index[mod1(rows[index], n_all)]
                if k > 0 && distance[k] == typemax(Int64)
                    distance[k] = level
                    push!(next, k)
                end
            end
        end
        frontier = next
    end
    return distance
end

"""
    level_counts(distance)

Number of nodes on every level of a distance vector, and the unreachable ones.
"""
function level_counts(distance::Vector{Int64})
    reached = filter(!=(typemax(Int64)), distance)
    return [count(==(l), reached) for l in 1:maximum(reached; init = 0)],
           length(distance) - length(reached)
end

"""
    wavefront_order(K0, nodes, dof, ordering)

Orders the condensed nodes for the cascade, so that consecutive chunks are compact and
the front stays narrow. The distance is counted in steps through the graph of `K0`.

- `:boundary`: from the farthest to the closest to the rest of the model, starting from
  every condensed node coupled to a node outside the condensed set. The chunks are
  layers parallel to the boundary, and the retained points join the front only in the
  last steps; but every layer is as large as the boundary, which for a thin region
  around a hole makes the front a whole ring.
- `:sweep`: Cuthill-McKee from a pseudo-peripheral node: the search starts at the
  condensed node farthest from the boundary, then once more from the node farthest from
  that one, and the nodes are taken in the order of that second search. The front is
  then a cross-section of the region moving through it, and the retained points join as
  it passes them.

Nodes without a path to the start come first; within a level the node ids decide.

# Arguments
- `K0::SparseMatrixCSC{Float64,Int64}`: Stiffness, only its sparsity is used
- `nodes::Vector{Int64}`: Condensed node ids, equal to their first degree of freedom
- `dof::Int64`: Degrees of freedom per node
- `ordering::Symbol`: `:boundary` or `:sweep`
# Returns
- `::Vector{Int64}`: Permutation of `eachindex(nodes)`
"""
function wavefront_order(K0::SparseMatrixCSC{Float64,Int64}, nodes::Vector{Int64},
                         dof::Int64, ordering::Symbol)
    n_all = size(K0, 1) ÷ dof
    local_index = zeros(Int64, n_all)
    local_index[nodes] .= eachindex(nodes)
    rows = rowvals(K0)
    values = nonzeros(K0)

    boundary = Int64[]
    for (j, node) in enumerate(nodes)
        coupled_outside = false
        for d in 1:dof, index in nzrange(K0, (d - 1) * n_all + node)
            if values[index] != 0.0 && local_index[mod1(rows[index], n_all)] == 0
                coupled_outside = true
                break
            end
        end
        coupled_outside && push!(boundary, j)
    end
    distance = graph_distance(K0, nodes, dof, local_index, boundary)

    if ordering == :boundary
        counts, unreached = level_counts(distance)
        @info "Craig-Bampton Cascade wavefront (boundary): $(length(counts)) levels from " *
              "the boundary inward, nodes per level $(counts), $unreached nodes without " *
              "a path to the boundary"
        return sortperm(collect(eachindex(nodes));
                        by = j -> (distance[j] == typemax(Int64) ? typemin(Int64) :
                                   -distance[j], nodes[j]))
    elseif ordering == :sweep
        reachable = findall(!=(typemax(Int64)), distance)
        start = isempty(reachable) ? 1 :
                reachable[argmax(distance[reachable])]
        first_pass = graph_distance(K0, nodes, dof, local_index, [start])
        reached = findall(!=(typemax(Int64)), first_pass)
        start = reached[argmax(first_pass[reached])]
        distance = graph_distance(K0, nodes, dof, local_index, [start])
        counts, unreached = level_counts(distance)
        @info "Craig-Bampton Cascade wavefront (sweep from node $(nodes[start])): " *
              "$(length(counts)) levels, nodes per level $(counts), $unreached nodes " *
              "not reached"
        return sortperm(collect(eachindex(nodes));
                        by = j -> (distance[j] == typemax(Int64) ? typemin(Int64) :
                                   distance[j], nodes[j]))
    end
    throw(ArgumentError("Unknown ordering $ordering, use :boundary or :sweep."))
end

"""
    split_subregions(c, dof, n_subregions, K0, ordering = :boundary)

Splits the condensed degrees of freedom into subregions of consecutive nodes in
wavefront order.

`c` holds the condensed degrees of freedom component by component, `[x of every node;
y of every node; ...]`, with the first component's index equal to the node id (see
`create_permutation`). The nodes are put in [`wavefront_order`](@ref) and cut into
`n_subregions` chunks of nearly equal size; every chunk takes all degrees of freedom of
its nodes, so a node is never split between two subregions. The result of the static
part is exact for any partition; the order only decides how wide the front gets, and
with it the memory.

# Arguments
- `c::Vector{Int64}`: Condensed degrees of freedom, component by component
- `dof::Int64`: Degrees of freedom per node
- `n_subregions::Int64`: Number of subregions, capped at the number of nodes
- `K0::SparseMatrixCSC{Float64,Int64}`: Stiffness, only its sparsity is used
- `ordering::Symbol`: Node order, see [`wavefront_order`](@ref)
# Returns
- `::Vector{Vector{Int64}}`: Degrees of freedom of every subregion
"""
function split_subregions(c::Vector{Int64}, dof::Int64, n_subregions::Int64,
                          K0::SparseMatrixCSC{Float64,Int64}, ordering::Symbol = :boundary)
    n_nodes = length(c) ÷ dof
    if n_nodes * dof != length(c)
        throw(ArgumentError("$(length(c)) condensed degrees of freedom are not a " *
                            "multiple of dof = $dof."))
    end
    n_subregions = clamp(n_subregions, 1, max(n_nodes, 1))
    order = wavefront_order(K0, c[1:n_nodes], dof, ordering)
    bounds = round.(Int64, range(0, n_nodes; length = n_subregions + 1))
    return [vec([c[(d-1)*n_nodes+j]
                 for j in order[(bounds[i] + 1):bounds[i + 1]], d in 1:dof])
            for i in 1:n_subregions]
end

"""
    plan_cascade(K0, subregions, active, n_modes)

Runs the front bookkeeping of the cascade on the sparsity of `K0` alone, before any
numbers are computed.

Which degrees of freedom join the front and which leave it depends only on which
entries of `K0` are nonzero, not on their values; every level condenses the previous
level's modes together with its subregion and adds up to `n_modes` new ones.
Playing it through once gives the new neighbours of every subregion, reused by the
actual cascade, and the largest sizes, so that every buffer can be allocated once.

# Arguments
- `K0::SparseMatrixCSC{Float64,Int64}`: Original stiffness
- `subregions::Vector{Vector{Int64}}`: Degrees of freedom of every subregion
- `active::BitVector`: Degrees of freedom taking part in the reduction
- `n_modes::Int64`: Modes computed per level
# Returns
- `neighbours::Vector{Vector{Int64}}`: Degrees of freedom joining the front with every
  subregion, sorted
- `max_ny::Int64`: Largest condensed set, subregion plus the previous level's modes
- `max_ns::Int64`: Largest coupling width, the front without the new modes
- `max_front::Int64`: Largest front including the new modes
"""
function plan_cascade(K0::SparseMatrixCSC{Float64,Int64},
                      subregions::Vector{Vector{Int64}}, active::BitVector,
                      n_modes::Int64)
    on_front = falses(length(active))
    eliminated = falses(length(active))
    n_physical = 0
    n_previous_modes = 0
    neighbours = Vector{Vector{Int64}}(undef, length(subregions))
    is_condensed = falses(length(active))
    for c in subregions
        is_condensed[c] .= true
    end
    max_ny = 0
    max_ns = 0
    max_front = 0
    for (i, c) in enumerate(subregions)
        candidates = union(nonzero_rows(K0[:, c]), nonzero_columns(K0[c, :]))
        eliminated[c] .= true
        n_physical -= count(view(on_front, c))
        on_front[c] .= false
        neighbours[i] = sort!([d
                               for d in candidates
                               if active[d] && !eliminated[d] && !on_front[d]])
        on_front[neighbours[i]] .= true
        n_physical += length(neighbours[i])
        ny = length(c) + n_previous_modes
        ns = n_physical
        nm = min(n_modes, ny)
        n_previous_modes = nm
        max_ny = max(max_ny, ny)
        max_ns = max(max_ns, ns)
        max_front = max(max_front, ns + nm)
        new_retained = count(d -> !is_condensed[d], neighbours[i])
        retained_on_front = count(d -> on_front[d] && !is_condensed[d],
                                  eachindex(on_front))
        # if i==19
        #     readline()
        # end
        @info "Craig-Bampton Cascade plan: subregion $i, nc=$(length(c)), ny=$ny, " *
              "coupling ns=$ns; new neighbours $(length(neighbours[i])) " *
              "($new_retained retained, $(length(neighbours[i]) - new_retained) " *
              "condensed later); retained on the front $retained_on_front, modes $nm"
    end
    return neighbours, max_ny, max_ns, max_front
end

"""
    Front

The part of the cascade's system that differs from the original matrices.

Condensing a subregion changes the stiffness and the mass only among the degrees of
freedom coupled to it, and replaces the modal coordinates of the previous level by its
own. Everything else keeps its
original entries. The cascade therefore never copies the original matrices: it keeps
them as they are and holds the changes as dense corrections on the front, the degrees of
freedom touched so far that are still in the system: the band of not yet condensed
points next to the condensed part, the retained points coupled to it, and the modal
coordinates of the current level. The effective matrices are the original ones plus these corrections.

The corrections live in the leading `length(labels)` rows and columns of `D_K` and
`D_M`, which are allocated once at the largest front size and reused by every step.

# Fields
- `labels::Vector{Int64}`: Global degree of freedom of every front entry, or `-m` for
  the m-th modal coordinate
- `position::Dict{Int64,Int64}`: Label => position in `labels`
- `D_K::Matrix{Float64}`: Stiffness correction, leading block in use
- `D_M::Matrix{Float64}`: Mass correction, leading block in use, symmetric
"""
mutable struct Front
    labels::Vector{Int64}
    position::Dict{Int64,Int64}
    D_K::Matrix{Float64}
    D_M::Matrix{Float64}
end

function Front(max_front::Int64)
    Front(Int64[], Dict{Int64,Int64}(), zeros(max_front, max_front),
          zeros(max_front, max_front))
end

"""
    Buffers

Work space of one cascade step, allocated once at the largest condensed set `max_ny`
(subregion plus the previous level's modes) and the largest coupling width `max_ns`;
every step works on leading views of it.

# Fields
- `Kcc::Matrix{Float64}`, `Mcc::Matrix{Float64}`, `C::Matrix{Float64}`: `max_ny` square
- `Kcs::Matrix{Float64}`, `T::Matrix{Float64}`: `max_ny x max_ns`; `Kcs` turns into `B`
- `Ksc::Matrix{Float64}`, `Msc::Matrix{Float64}`: `max_ns x max_ny`
"""
struct Buffers
    Kcc::Matrix{Float64}
    Mcc::Matrix{Float64}
    C::Matrix{Float64}
    Kcs::Matrix{Float64}
    T::Matrix{Float64}
    Ksc::Matrix{Float64}
    Msc::Matrix{Float64}
end

function Buffers(max_ny::Int64, max_ns::Int64)
    Buffers(zeros(max_ny, max_ny), zeros(max_ny, max_ny), zeros(max_ny, max_ny),
            zeros(max_ny, max_ns), zeros(max_ny, max_ns), zeros(max_ns, max_ny),
            zeros(max_ns, max_ny))
end

"""
    compact!(D, keep)

Moves rows and columns `keep` of `D`, ascending, to its leading block in place.

Every entry moves up and to the left or stays, so walking the target in column major
order never overwrites an entry that is still to be read.
"""
function compact!(D::Matrix{Float64}, keep::Vector{Int64})
    @inbounds for (j_new, j_old) in enumerate(keep), (i_new, i_old) in enumerate(keep)
        D[i_new, j_new] = D[i_old, j_old]
    end
    return D
end

"""
    condense_subregion!(front, buffers, K0, M0, c, neighbours, eliminated, n_modes,
                        first_mode, max_frequency)

One level of the cascade: condenses the subregion together with the previous level's
modal coordinates, keeping up to `n_modes` new fixed-interface modes, and updates the
front in place.

The condensed set is `y = [c; previous modes]`. The effective matrices are `K0`,
`diag(M0)` plus the front corrections; the previous modes have no part in `K0` and `M0`
and live in the front alone. With `s` the degrees of freedom coupled to `y` -- the
physical front without the subregion, plus its new neighbours -- `y` is represented by
its static recovery `-B u_s`, `B = Kyy^-1 Kys`, plus the fixed-interface modes `X` of
`Kyy X = Myy X W` (written with `c` for `y` below):

    Kss' = Kss - Ksc B
    Ksm  = Ksc X - B' Kcc X       Kms = X' Kcs - X' Kcc B       Kmm = X' Kcc X
    Mss' = Mss - Msc B - B' Msc' + B' Mcc B
    Msm  = Msc X - B' Mcc X                                     Mmm = X' Mcc X

The mass corrections are symmetric, so `Mcs = Msc'`. The original part of `Kss` and
`Mss` stays in `K0` and `M0`; only the changes go into the front. The stiffness coupling
blocks are formed explicitly, since the peridynamic stiffness need not be symmetric (see
`CraigBampton.reduce_matrices`).

The modes come from the standard problem `L^-1 Kcc L^-T V = V W` with `Mcc = L L'`,
of which only the `n_modes` lowest eigenpairs are computed; `X = L^-T V` is mass
normalised. With `max_frequency` set, those above it are dropped, the same cutoff on
every level.

Nothing of the size of the subregion or the front is allocated: every block is a view
of `buffers` or of the front.

# Arguments
- `front::Front`: Front corrections, updated in place
- `buffers::Buffers`: Work space
- `K0::SparseMatrixCSC{Float64,Int64}`: Original stiffness, not modified
- `M0::Vector{Float64}`: Original lumped mass, not modified
- `c::Vector{Int64}`: Global degrees of freedom of the subregion
- `neighbours::Vector{Int64}`: Degrees of freedom joining the front, from
  [`plan_cascade`](@ref)
- `eliminated::BitVector`: Degrees of freedom condensed so far, updated in place
- `n_modes::Int64`: Number of fixed-interface modes to keep
- `first_mode::Int64`: Number of the first new modal coordinate
- `max_frequency::Union{Nothing,Float64}`: Highest mode frequency in Hz to keep
# Returns
- `w::Vector{Float64}`: Eigenvalues of the new modes, ascending
"""
function condense_subregion!(front::Front, buffers::Buffers,
                             K0::SparseMatrixCSC{Float64,Int64}, M0::Vector{Float64},
                             c::Vector{Int64}, neighbours::Vector{Int64},
                             eliminated::BitVector, n_modes::Int64, first_mode::Int64,
                             max_frequency::Union{Nothing,Float64})
    nc = length(c)
    @timeit "front bookkeeping" begin
        # condensed set y: the subregion, then the previous level's modes
        c_local = Int64[]
        c_front = Int64[]
        for (i, d) in enumerate(c)
            p = get(front.position, d, 0)
            if p > 0
                push!(c_local, i)
                push!(c_front, p)
            end
        end
        previous_modes = findall(<(0), front.labels)
        ny = nc + length(previous_modes)
        append!(c_local, (nc + 1):ny)
        append!(c_front, previous_modes)
        in_c = Set(c)
        kept_front = [p
                      for (p, label) in enumerate(front.labels)
                      if label > 0 && !(label in in_c)]
        nf = length(kept_front)
        s = vcat(front.labels[kept_front], neighbours)
        ns = length(s)
        s_physical = findall(>(0), s)
        nm = min(n_modes, ny)
    end
    @info "Craig-Bampton Cascade: subregion nc=$nc, ny=$ny, coupling ns=$ns ($nf " *
          "carried, $(length(neighbours)) new), front after $(ns + nm)"

    pc = 1:nc
    cc = 1:ny
    sc = 1:ns
    Kcc = view(buffers.Kcc, cc, cc)
    Mcc = view(buffers.Mcc, cc, cc)
    C = view(buffers.C, cc, cc)
    Kcs = view(buffers.Kcs, cc, sc)
    T = view(buffers.T, cc, sc)
    Ksc = view(buffers.Ksc, sc, cc)
    Msc = view(buffers.Msc, sc, cc)

    @timeit "fill blocks" begin
        fill!(Kcc, 0.0)
        fill!(Mcc, 0.0)
        fill!(Kcs, 0.0)
        fill!(Ksc, 0.0)
        fill!(Msc, 0.0)
        add_sparse!(Kcc, K0, c, c, pc, pc)
        add_sparse!(Kcs, K0, c, s[s_physical], pc, s_physical)
        add_sparse!(Ksc, K0, s[s_physical], c, s_physical, pc)
        @inbounds for i in pc
            Mcc[i, i] = M0[c[i]]
        end
        if !isempty(c_front)
            # the front without the subregion is the first nf entries of s
            add_block!(Kcc, c_local, c_local, front.D_K, c_front, c_front)
            add_block!(Mcc, c_local, c_local, front.D_M, c_front, c_front)
            add_block!(Kcs, c_local, 1:nf, front.D_K, c_front, kept_front)
            add_block!(Ksc, 1:nf, c_local, front.D_K, kept_front, c_front)
            add_block!(Msc, 1:nf, c_local, front.D_M, kept_front, c_front)
        end
    end

    @timeit "fixed interface modes" begin
        if nm > 0
            copyto!(C, Kcc)
            L = cholesky!(Symmetric(Mcc, :L)).L       # Mcc now holds its factor
            ldiv!(L, C)
            rdiv!(C, transpose(L))
            modes = eigen!(Symmetric(C, :L), 1:nm)
            w = modes.values
            X = modes.vectors
            ldiv!(transpose(L), X)
            if !isnothing(max_frequency)
                n_below = count(<=(max_frequency), sqrt.(max.(w, 0.0)) ./ (2 * pi))
                w = w[1:n_below]
                X = X[:, 1:n_below]
            end
        else
            L = cholesky!(Symmetric(Mcc, :L)).L
            w = Float64[]
            X = zeros(ny, 0)
        end
        nm = length(w)
        n_new = ns + nm
    end
    @timeit "modal products" begin
        Kcc_X = Kcc * X
        Xt_Kcc = transpose(X) * Kcc
        Xt_Kcs = transpose(X) * Kcs
        Mcc_X = L * (transpose(L) * X)
    end

    # Kcc is overwritten by its factorization, Kcs by the recovery modes B.
    @timeit "factorize Kcc" factorization=factorize_condensed!(Kcc)
    @timeit "recovery modes" B=ldiv!(factorization, Kcs)
    @timeit "Mcc B" begin
        copyto!(T, B)
        lmul!(transpose(L), T)
        lmul!(L, T)
    end

    @timeit "compact front" begin
        compact!(front.D_K, kept_front)
        compact!(front.D_M, kept_front)
        front.D_K[(nf + 1):n_new, 1:n_new] .= 0.0
        front.D_K[1:nf, (nf + 1):n_new] .= 0.0
        front.D_M[(nf + 1):n_new, 1:n_new] .= 0.0
        front.D_M[1:nf, (nf + 1):n_new] .= 0.0
    end

    se = 1:ns
    me = (ns + 1):n_new
    D_K = front.D_K
    D_M = front.D_M
    @timeit "front stiffness" begin
        @views mul!(D_K[se, se], Ksc, B, -1.0, 1.0)
        if nm > 0
            @views mul!(D_K[se, me], Ksc, X)
            @views mul!(D_K[se, me], transpose(B), Kcc_X, -1.0, 1.0)
            @views D_K[me, se] .= Xt_Kcs
            @views mul!(D_K[me, se], Xt_Kcc, B, -1.0, 1.0)
            @views mul!(D_K[me, me], Xt_Kcc, X)
        end
    end
    @timeit "front mass" begin
        @views mul!(D_M[se, se], transpose(B), T, 1.0, 1.0)
        @views mul!(D_M[se, se], Msc, B, -1.0, 1.0)
        @views mul!(D_M[se, se], transpose(B), transpose(Msc), -1.0, 1.0)
        if nm > 0
            @views mul!(D_M[se, me], Msc, X)
            @views mul!(D_M[se, me], transpose(B), Mcc_X, -1.0, 1.0)
            @views D_M[me, se] .= transpose(D_M[se, me])
            @views mul!(D_M[me, me], transpose(X), Mcc_X)
        end
    end

    @timeit "front update" begin
        front.labels = vcat(s, .-(first_mode:(first_mode + nm - 1)))
        front.position = Dict(label => p for (p, label) in enumerate(front.labels))
        eliminated[c] .= true
    end
    return w
end

"""
    add_block!(target, target_rows, target_columns, source, rows, columns)

Adds `source[rows, columns]` into `target[target_rows, target_columns]` entry by entry.
A broadcast over the two indexed views would copy the source first whenever it cannot
rule out that both share memory.
"""
function add_block!(target::AbstractMatrix{Float64}, target_rows::AbstractVector{Int64},
                    target_columns::AbstractVector{Int64}, source::Matrix{Float64},
                    rows::AbstractVector{Int64}, columns::AbstractVector{Int64})
    @inbounds for (tj, j) in zip(target_columns, columns), (ti, i) in zip(target_rows,
                                                                          rows)
        target[ti, tj] += source[i, j]
    end
    return target
end

"""
    add_sparse!(target, A, rows, columns, target_rows, target_columns)

Adds `A[rows, columns]` into `target[target_rows, target_columns]` without the
temporary the slice would allocate. `rows` must be sorted.
"""
function add_sparse!(target::AbstractMatrix{Float64}, A::SparseMatrixCSC{Float64,Int64},
                     rows::AbstractVector{Int64}, columns::AbstractVector{Int64},
                     target_rows::AbstractVector{Int64},
                     target_columns::AbstractVector{Int64})
    A_rows = rowvals(A)
    A_values = nonzeros(A)
    order = sortperm(rows)
    sorted_rows = rows[order]
    @inbounds for (k, column) in enumerate(columns)
        tc = target_columns[k]
        for index in nzrange(A, column)
            row = A_rows[index]
            i = searchsortedfirst(sorted_rows, row)
            if i <= length(sorted_rows) && sorted_rows[i] == row
                target[target_rows[order[i]], tc] += A_values[index]
            end
        end
    end
    return target
end

"""
    merge_front_block(base, D, positions, n)

The reduced matrix `base + D`, assembled column by column straight into its final
compressed sparse column form.

`base` is the original part, `n0 x n0` with `n0 <= n`; `D` is the dense front block whose
rows and columns go to `positions`. Where both have an entry they are summed. Neither a
sparse copy of `D` nor a second matrix for the sum is built: the result is the only
allocation of its size.

# Arguments
- `base::SparseMatrixCSC{Float64,Int64}`: Original part, the leading `n0 x n0` block
- `D::AbstractMatrix{Float64}`: Dense front block
- `positions::Vector{Int64}`: Row and column of every entry of `D` in the result
- `n::Int64`: Size of the result
# Returns
- `::SparseMatrixCSC{Float64,Int64}`: The reduced matrix
"""
function merge_front_block(base::SparseMatrixCSC{Float64,Int64}, D::AbstractMatrix{Float64},
                           positions::Vector{Int64}, n::Int64)
    order = sortperm(positions)
    sorted = positions[order]
    m = length(sorted)
    n0 = size(base, 2)
    base_rows = rowvals(base)
    base_values = nonzeros(base)
    # column of the result => index into sorted, 0 if the column has no block entries
    block_column = zeros(Int64, n)
    block_column[sorted] .= eachindex(sorted)
    in_block = falses(n)
    in_block[sorted] .= true

    colptr = Vector{Int64}(undef, n + 1)
    colptr[1] = 1
    @inbounds for column in 1:n
        count_column = block_column[column] > 0 ? m : 0
        if column <= n0
            for index in nzrange(base, column)
                (block_column[column] > 0 && in_block[base_rows[index]]) && continue
                count_column += 1
            end
        end
        colptr[column + 1] = colptr[column] + count_column
    end

    rowval = Vector{Int64}(undef, colptr[n + 1] - 1)
    nzval = Vector{Float64}(undef, colptr[n + 1] - 1)
    @inbounds for column in 1:n
        k = colptr[column]
        jj = block_column[column]
        base_range = column <= n0 ? nzrange(base, column) : (1:0)
        if jj == 0
            for index in base_range
                rowval[k] = base_rows[index]
                nzval[k] = base_values[index]
                k += 1
            end
            continue
        end
        # merge the sorted base rows with the sorted block rows
        j = order[jj]
        b = first(base_range)
        b_end = last(base_range)
        for ii in 1:m
            row = sorted[ii]
            while b <= b_end && base_rows[b] < row
                rowval[k] = base_rows[b]
                nzval[k] = base_values[b]
                k += 1
                b += 1
            end
            value = D[order[ii], j]
            if b <= b_end && base_rows[b] == row
                value += base_values[b]
                b += 1
            end
            rowval[k] = row
            nzval[k] = value
            k += 1
        end
        while b <= b_end
            rowval[k] = base_rows[b]
            nzval[k] = base_values[b]
            k += 1
            b += 1
        end
    end
    return SparseMatrixCSC(n, n, colptr, rowval, nzval)
end

"""
    reduce_matrices(K, M_diag, r, c, n_modes = 0; check_symmetry_sample = 2000, dof = 1,
                    n_subregions = 20, max_frequency = nothing, kwargs...)

Craig-Bampton reduction done as a cascade over subregions of the condensed region,
instead of one factorization of the whole region.

The condensed points are split into `n_subregions` layers in wavefront order, from the
farthest to the closest to the rest of the model ([`split_subregions`](@ref)) and condensed one after the other
([`condense_subregion!`](@ref)). The original matrices are never copied or rebuilt: the
changes the condensation makes are held as dense corrections on the front
([`Front`](@ref)). A first pass over the sparsity alone ([`plan_cascade`](@ref)) gives
the largest front and subregion, so that the front and all work space are allocated once
and reused by every step.

The static part is exact: condensing one subregion after the other gives the same Schur
complement as condensing all of them at once, for any partition. `n_modes = 0` is
therefore Guyan condensation of the whole region.

The modes are passed on from level to level: every level condenses the previous
level's modal coordinates together with its subregion and computes up to `n_modes` new
fixed-interface modes of that combined set, so they describe everything condensed so
far. Only the modes of the last level remain. With `max_frequency` set, every level
keeps only the modes up to that frequency, the same cutoff throughout. The truncation
errors of the levels accumulate.

# Arguments
- `K::AbstractMatrix`: Stiffness matrix
- `M_diag::AbstractVector`: Lumped mass matrix as a vector, one entry per degree of freedom
- `r::AbstractVector{<:Integer}`: Indices of the retained degrees of freedom
- `c::AbstractVector{<:Integer}`: Indices of the condensed degrees of freedom, component by
  component (see [`split_subregions`](@ref))
- `n_modes::Integer`: Number of modes per level, and so in the result; 0 is Guyan
  condensation
# Keywords
- `check_symmetry_sample::Int64`: Entries drawn for the symmetry check, 0 disables it
- `dof::Int64`: Degrees of freedom per node
- `n_subregions::Int64`: Number of subregions the condensed region is split into
- `max_frequency::Union{Nothing,Float64}`: Highest mode frequency in Hz to keep
# Returns
- `K_reduced::SparseMatrixCSC`: Reduced stiffness, `r` in the given order followed by the
  modes of the last level in ascending frequency
- `M_reduced::SparseMatrixCSC`: Reduced mass, same layout
"""
function reduce_matrices(K::AbstractMatrix,
                         M_diag::AbstractVector,
                         r::AbstractVector{<:Integer},
                         c::AbstractVector{<:Integer},
                         n_modes::Integer = 0;
                         check_symmetry_sample::Int64 = 2000,
                         dof::Int64 = 1,
                         n_subregions::Int64 = 20,
                         max_frequency::Union{Nothing,Float64} = nothing,
                         kwargs...)
    n_modes = Int64(n_modes)
    n_modes < 0 && throw(ArgumentError("n_modes = $n_modes, must be >= 0."))
    r = collect(Int64, r)
    c = collect(Int64, c)
    if !isempty(intersect(r, c))
        throw(ArgumentError("Retained and condensed index sets overlap."))
    end
    K0 = K isa SparseMatrixCSC{Float64,Int64} ? K :
         SparseMatrixCSC{Float64,Int64}(sparse(K))
    M0 = collect(Float64, M_diag)

    if check_symmetry_sample > 0
        @timeit "cascade symmetry check" check_symmetry(K0; samples = check_symmetry_sample)
    end

    active = falses(size(K0, 1))
    active[r] .= true
    active[c] .= true

    # Thinner layers take band points off the front. If the buffers do not fit, the
    # number of subregions is doubled and the plan repeated, as long as the buffers still
    # shrink noticeably and the subregions keep at least 10 nodes; the retained boundary
    # itself is a lower bound no partition can undercut.
    n_nodes = length(c) ÷ dof
    previous_mib = Inf
    local subregions, neighbours, max_ny, max_ns, max_front, buffer_mib, ordering
    while true
        # both orders are planned, the one with the smaller buffers is used
        local best = nothing
        for ordering in (:boundary, :sweep)
            @timeit "split subregions" candidate=split_subregions(c, dof, n_subregions,
                                                                  K0, ordering)
            @timeit "plan cascade" plan=plan_cascade(K0, candidate, active, n_modes)
            candidate_ny, candidate_ns, candidate_front = plan[2], plan[3], plan[4]
            candidate_mib = 8 * (2 * candidate_front^2 + 3 * candidate_ny^2 +
                             4 * candidate_ny * candidate_ns) / 2^20
            @info "Craig-Bampton Cascade: ordering $ordering, front $candidate_front, " *
                  "buffers $(round(candidate_mib; digits = 1)) MiB"
            if isnothing(best) || candidate_mib < best[2]
                best = ((ordering, candidate, plan), candidate_mib)
            end
        end
        (ordering, subregions, (neighbours, max_ny, max_ns, max_front)), buffer_mib = best
        # Garbage from the setup and from earlier plans counts as used memory until
        # collected.
        @timeit "garbage collection" GC.gc()
        available_mib = Sys.free_memory() / 2^20
        @info "Craig-Bampton Cascade: $(length(r)) retained, $(length(c)) condensed " *
              "degrees of freedom in $(length(subregions)) subregions ($ordering order), " *
              "$n_modes modes; largest condensed set $max_ny, coupling $max_ns, front " *
              "$max_front; buffers " *
              "$(round(buffer_mib; digits = 1)) MiB, $(round(available_mib; digits = 1)) " *
              "MiB free"
        buffer_mib <= available_mib && break

        if n_nodes ÷ (2 * length(subregions)) < 10 || buffer_mib > 0.95 * previous_mib
            error("Craig-Bampton Cascade: the buffers need " *
                  "$(round(buffer_mib; digits = 1)) MiB, but only " *
                  "$(round(available_mib; digits = 1)) MiB are free, and more subregions " *
                  "no longer shrink them; the front reaches $max_front degrees of freedom. See the " *
                  "plans above: the retained points on the front are the boundary of the " *
                  "condensed region, which no partition can make smaller.")
        end
        previous_mib = buffer_mib
        n_subregions = 2 * length(subregions)
        @info "Craig-Bampton Cascade: buffers do not fit, retrying with $n_subregions " *
              "subregions"
    end
    @timeit "allocate buffers" begin
        front = Front(max_front)
        buffers = Buffers(max_ny, max_ns)
    end
    eliminated = falses(size(K0, 1))
    frequencies = Float64[]

    for (i, subregion) in enumerate(subregions)
        w = @timeit "cascade subregion" condense_subregion!(front, buffers, K0, M0,
                                                            subregion, neighbours[i],
                                                            eliminated, n_modes,
                                                            length(frequencies) + 1,
                                                            max_frequency)
        append!(frequencies, sqrt.(max.(w, 0.0)) ./ (2 * pi))
        @info "Craig-Bampton Cascade: subregion $i of $(length(subregions)), " *
              "$(length(subregion)) condensed, $(length(w)) modes, front " *
              "$(length(front.labels)), peak RSS " *
              "$(round(Sys.maxrss() / 2^20; digits = 1)) MiB"
        @timeit "garbage collection" GC.gc(false)
    end
    buffers = nothing
    @timeit "garbage collection" GC.gc()

    # the modes of the last level remain; they come in ascending frequency
    selected = findall(<(0), front.labels)
    if !isnothing(max_frequency) && n_modes > 0 && length(selected) == n_modes
        @warn "Craig-Bampton Cascade: all $(length(selected)) modes of the last level lie " *
              "below the maximum frequency of $(max_frequency) Hz; raise " *
              "'Number of Modes' to find more."
    end

    nr = length(r)
    n_kept = length(selected)
    n_block = count(>(0), front.labels) + n_kept
    # the front block goes into both reduced matrices as sparse entries, 16 bytes each
    result_mib = 2 * 16 * n_block^2 / 2^20
    available_mib = Sys.free_memory() / 2^20
    @info "Craig-Bampton Cascade: dense boundary block $n_block x $n_block, assembling " *
          "the reduced matrices needs about $(round(result_mib; digits = 1)) MiB, " *
          "$(round(available_mib; digits = 1)) MiB free"
    if result_mib > available_mib
        error("Craig-Bampton Cascade: the reduced matrices need about " *
              "$(round(result_mib; digits = 1)) MiB, but only " *
              "$(round(available_mib; digits = 1)) MiB are free. Their boundary block is " *
              "dense, $n_block x $n_block: the retained degrees of freedom coupled to " *
              "the condensed region plus the modes.")
    end
    @timeit "final reduced matrices" begin
        r_index = Dict(g => i for (i, g) in enumerate(r))
        target = zeros(Int64, length(front.labels))
        for (p, label) in enumerate(front.labels)
            label > 0 && (target[p] = r_index[label])
        end
        for (k, p) in enumerate(selected)
            target[p] = nr + k
        end
        kept = findall(>(0), target)
        n_total = nr + n_kept
        K_reduced = merge_front_block(K0[r, r], view(front.D_K, kept, kept),
                                      target[kept], n_total)
        M_reduced = merge_front_block(spdiagm(M0[r]), view(front.D_M, kept, kept),
                                      target[kept], n_total)
        dropzeros!(K_reduced)
        dropzeros!(M_reduced)
    end

    if n_kept > 0
        kept_frequencies = frequencies[.-front.labels[selected]]
        @info "Craig-Bampton Cascade: $n_kept modes of the last level kept"
        @info "Craig-Bampton fixed-interface frequencies: " *
              "$(round(kept_frequencies[1]; sigdigits = 4)) Hz (mode 1) to " *
              "$(round(kept_frequencies[end]; sigdigits = 4)) Hz (mode $n_kept)"
        @info "Craig-Bampton fixed-interface frequency list [Hz]: " *
              join(round.(kept_frequencies; sigdigits = 6), ", ")
    end
    report_modal_coupling(K_reduced, M_reduced, nr, n_kept)

    return K_reduced, M_reduced
end

end
