# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Model_reduction
using TimerOutputs: @timeit
using ...Data_Manager
using SparseArrays
using Serialization
using ....ModuleLoader: find_module_files, create_module_specifics
using ...Correspondence_matrix_based: build_mass_matrix, init_model, init_matrix,
                                      rebuild_matrix!
using ...Helpers: find_active_nodes, create_permutation
global module_list = find_module_files(@__DIR__, "model_reduction_name")
for mod in module_list
    include(mod["File"])
end

export ReducedState
export n_physical, n_modal, physical_part, modal_part
export pull_from_nodes!, add_from_nodes!, push_to_nodes!
export setup_reduced_state

"""
    ReducedState

State vectors of the time integration, held flat rather than as node fields.

The vector is ordered `[u_retained; eta]`: `n_phys` physical degrees of freedom of the
retained nodes, followed by `n_modal` modal coordinates. A modal coordinate belongs to no
node and therefore cannot live in a `nnodes x dof` field, which is what this container
is for.

`n_modal` is zero for a static reduction (Guyan) and for a run without any reduction, in
which case the state is nothing but the retained degrees of freedom and the integration
is the same as before. The unreduced case is the one where every node is retained.

The physical part is ordered component by component, the same way `vec` flattens an
`n x dof` node field: first the first component of every retained node, then the second,
and so on. Retained node `k` therefore sits at `(d - 1) * n_retained + k` for component
`d`. This has to match the ordering the reduced matrices were built in.

# Fields
- `n_phys::Int64`: Number of physical degrees of freedom
- `n_modal::Int64`: Number of modal coordinates
- `q::Vector{Float64}`: Displacements
- `q_dot::Vector{Float64}`: Velocities
- `q_ddot::Vector{Float64}`: Accelerations
- `f::Vector{Float64}`: Force scratch vector, reused every step
"""
struct ReducedState
    n_phys::Int64
    n_modal::Int64
    q::Vector{Float64}
    q_dot::Vector{Float64}
    q_ddot::Vector{Float64}
    f::Vector{Float64}
end

"""
    ReducedState(n_phys::Int64, n_total::Int64)

Allocates the state for a system of `n_total` degrees of freedom, the first `n_phys` of
them physical.

`n_total` comes from `size(K, 1)`, so the number of modal coordinates never has to be
passed around: it is whatever the reduction scheme added.

# Arguments
- `n_phys::Int64`: Number of physical degrees of freedom
- `n_total::Int64`: Size of the system
# Returns
- `::ReducedState`: Zero initialised state
"""
function ReducedState(n_phys::Int64, n_total::Int64)
    if n_total < n_phys
        throw(ArgumentError("System has $n_total degrees of freedom, fewer than the " *
                            "$n_phys physical ones."))
    end
    return ReducedState(n_phys, n_total - n_phys, zeros(n_total), zeros(n_total),
                        zeros(n_total), zeros(n_total))
end

"""
    n_physical(state)

Number of physical degrees of freedom.
"""
n_physical(state::ReducedState) = state.n_phys

"""
    n_modal(state)

Number of modal coordinates; zero without a modal reduction.
"""
n_modal(state::ReducedState) = state.n_modal

"""
    physical_part(vector, state)

View on the entries that belong to retained nodes.

Everything with a physical meaning — boundary conditions, external loads, output — must
only touch this part.
"""
physical_part(v::AbstractVector, state::ReducedState) = @view v[1:state.n_phys]

"""
    modal_part(vector, state)

View on the modal coordinates; empty without a modal reduction.
"""
modal_part(v::AbstractVector, state::ReducedState) = @view v[(state.n_phys + 1):end]

"""
    pull_from_nodes!(vector, field, retained_nodes, state)

Copies the retained rows of a node field into the physical part of a state vector.

The modal part is left untouched, because a node field says nothing about it.

# Arguments
- `vector::AbstractVector{Float64}`: Target, length `n_phys + n_modal`
- `field::AbstractMatrix{Float64}`: Node field, `nnodes x dof`
- `retained_nodes::AbstractVector{Int64}`: Retained node indices
- `state::ReducedState`: The state, for the split position
# Returns
- `vector`: The unchanged reference
"""
function pull_from_nodes!(vector::AbstractVector{Float64},
                          field::AbstractMatrix{Float64},
                          retained_nodes::AbstractVector{Int64},
                          state::ReducedState)
    dof = size(field, 2)
    n_retained = length(retained_nodes)
    # Component major, matching vec(field[retained_nodes, :]).
    @inbounds for d in 1:dof
        offset = (d - 1) * n_retained
        for (k, node) in enumerate(retained_nodes)
            vector[offset+k] = field[node, d]
        end
    end
    return vector
end

"""
    add_from_nodes!(vector, field, retained_nodes, state)

Adds the retained rows of a node field onto the physical part of a state vector.

Same as [`pull_from_nodes!`](@ref), but accumulates instead of overwriting -- for adding
a node field's contribution onto a vector that already holds something else, such as
`state.f` already holding the stiffness-matrix force before the material point
contribution is added onto it. The modal part is left untouched, for the same reason as
`pull_from_nodes!`.

# Arguments
- `vector::AbstractVector{Float64}`: Target, length `n_phys + n_modal`, mutated in place
- `field::AbstractMatrix{Float64}`: Node field, `nnodes x dof`
- `retained_nodes::AbstractVector{Int64}`: Retained node indices
- `state::ReducedState`: The state, for the split position
# Returns
- `vector`: The same vector, mutated
"""
function add_from_nodes!(vector::AbstractVector{Float64},
                         field::AbstractMatrix{Float64},
                         retained_nodes::AbstractVector{Int64},
                         state::ReducedState)
    dof = size(field, 2)
    n_retained = length(retained_nodes)
    # Component major, matching vec(field[retained_nodes, :]).
    @inbounds for d in 1:dof
        offset = (d - 1) * n_retained
        for (k, node) in enumerate(retained_nodes)
            vector[offset+k] += field[node, d]
        end
    end
    return vector
end

"""
    push_to_nodes!(field, vector, retained_nodes, state)

Copies the physical part of a state vector back into the retained rows of a node field.

Rows of nodes that are not retained are left alone, and the modal part is dropped: it has
no node to be written to.

# Arguments
- `field::AbstractMatrix{Float64}`: Node field, `nnodes x dof`
- `vector::AbstractVector{Float64}`: Source, length `n_phys + n_modal`
- `retained_nodes::AbstractVector{Int64}`: Retained node indices
- `state::ReducedState`: The state, for the split position
# Returns
- `field`: The unchanged reference
"""
function push_to_nodes!(field::AbstractMatrix{Float64},
                        vector::AbstractVector{Float64},
                        retained_nodes::AbstractVector{Int64},
                        state::ReducedState)
    dof = size(field, 2)
    n_retained = length(retained_nodes)
    # Component major, matching vec(field[retained_nodes, :]).
    @inbounds for d in 1:dof
        offset = (d - 1) * n_retained
        for (k, node) in enumerate(retained_nodes)
            field[node, d] = vector[offset+k]
        end
    end
    return field
end

"""
    setup_reduced_state(model_reduction, K)

Retained nodes and state vectors for the time loop.

Without a reduction every node is retained and there are no modal coordinates, so the
state is the full system written as a flat vector and the time loop needs no branch on
whether a reduction is active.

# Arguments
- `model_reduction`: The `"Model Reduction"` solver option, `false` when disabled
- `K::AbstractMatrix`: The stiffness matrix, reduced or not
# Returns
- `retained_nodes::Vector{Int64}`: Nodes carrying physical degrees of freedom
- `state::ReducedState`: Zero initialised state
"""
function setup_reduced_state(model_reduction, K::AbstractMatrix)
    dof = Data_Manager.get_dof()

    retained_nodes = model_reduction == false ?
                     collect(1:Data_Manager.get_nnodes()) :
                     Data_Manager.get_reduced_model_retained()

    n_phys = length(retained_nodes) * dof
    n_total = size(K, 1)

    if n_total < n_phys
        throw(ArgumentError("Stiffness matrix has $n_total degrees of freedom for " *
                            "$n_phys retained degrees of freedom. Retained node list " *
                            "and matrix disagree."))
    end

    state = ReducedState(n_phys, n_total)
    if n_modal(state) > 0
        @info "Reduced system: $n_phys physical and $(n_modal(state)) modal degrees of freedom"
    end

    return retained_nodes, state
end

"""
    parse_reduction_blocks(model_param)

Block IDs to condense away, read from the `"Reduction Blocks"` input deck entry.

Accepts a single integer or a string of block IDs separated by commas/whitespace -- the
two forms the YAML parser can hand back. Logs and returns `nothing` if the entry is
absent or of an unsupported type, so the caller only has to check for that.

# Arguments
- `model_param::AbstractDict`: The `"Model Reduction"` solver parameters
# Returns
- `Union{Nothing,Vector{Int64}}`: The block IDs, or `nothing` if none were usable
"""
function parse_reduction_blocks(model_param::AbstractDict)
    reduction_blocks = get(model_param, "Reduction Blocks", nothing)
    if isnothing(reduction_blocks)
        @warn "No reduction blocks defined for model reduction. If you want to use a reduced model please define 'Reduction Blocks' in the yaml input deck."
        return nothing
    end
    if reduction_blocks isa Float64
        @error "Type Float is not supported for Reduction Blocks"
        return nothing
    end
    reduction_blocks isa Int64 && return [reduction_blocks]
    return parse.(Int64, filter(!isempty, split(reduction_blocks, r"[,\s]")))
end

"""
    partition_nodes(block_nodes, reduction_blocks, material_point_region)

Splits every node into the sets the reduction needs.

Retained nodes are the material point nodes plus every one of their bonded neighbours --
peridynamics' nonlocal (horizon-based) interactions mean that boundary layer has to stay
physical, even where it geometrically belongs to a reduction block. Coupling nodes are
exactly that layer, `retained \\ material point`: retained nodes with no separate force
computation of their own, relying entirely on the reduced operator.

`material_point_region = false` empties the material point set, folding every retained
node into the coupling layer instead.

# Arguments
- `block_nodes::AbstractDict{Int64,Vector{Int64}}`: Nodes per block
- `reduction_blocks::Vector{Int64}`: Block IDs to condense away
- `material_point_region::Bool`: Whether the non-reduction blocks keep their own force computation
# Returns
- `retained_nodes::Vector{Int64}`, `condensed_nodes::Vector{Int64}`,
  `pd_nodes::Vector{Int64}`, `coupling_nodes::Vector{Int64}`, all sorted
"""
function partition_nodes(block_nodes::AbstractDict{Int64,Vector{Int64}},
                         reduction_blocks::Vector{Int64}, material_point_region::Bool)
    nlist = Data_Manager.get_nlist()
    nnodes = Data_Manager.get_nnodes()

    full_blocks = setdiff(collect(keys(block_nodes)), reduction_blocks)
    pd_nodes = Int64[]
    for block in full_blocks
        append!(pd_nodes, block_nodes[block])
    end

    retained_nodes = Int64[]
    for node in pd_nodes
        append!(retained_nodes, nlist[node])
    end
    append!(retained_nodes, pd_nodes)
    retained_nodes = sort(unique(retained_nodes))

    material_point_region || empty!(pd_nodes)
    sort!(pd_nodes)

    condensed_nodes = sort!(setdiff(collect(1:nnodes), retained_nodes))
    coupling_nodes = setdiff(retained_nodes, pd_nodes)

    return retained_nodes, condensed_nodes, pd_nodes, coupling_nodes
end

"""
    mark_coupling_field!(retained_nodes, condensed_nodes, pd_nodes, coupling_nodes)

Fills the `"Coupling Nodes"` field for visualisation and debugging; it has no effect on
the reduction itself.
"""
function mark_coupling_field!(retained_nodes::Vector{Int64},
                              condensed_nodes::Vector{Int64},
                              pd_nodes::Vector{Int64}, coupling_nodes::Vector{Int64})
    cn = Data_Manager.create_constant_node_scalar_field("Coupling Nodes", Int64)
    cn[retained_nodes] .= 3
    cn[condensed_nodes] .= 6
    cn[pd_nodes] .+= 1  # added to be sure that all points are handled
    cn[coupling_nodes] .+= 3  # added to be sure that all points are handled
    return nothing
end

"""
    expand_density_per_dof(density, dof)

Repeats each node's density across its `dof` degrees of freedom, the lumped mass
convention `Matrix_Verlet.jl` also uses for the unreduced system.
"""
function expand_density_per_dof(density, dof::Int64)
    density_mass = zeros(Float64, length(density) * dof)
    for iID in eachindex(density_mass)
        density_mass[iID] = density[Int(ceil(iID / dof))]
    end
    return density_mass
end

function init_reduce_model(model_param::AbstractDict,
                           block_nodes::AbstractDict{Int64,Vector{Int64}},
                           density)
    reduction_blocks = parse_reduction_blocks(model_param)
    isnothing(reduction_blocks) && return

    @info "Model Reduction Type: $(model_param["Type"])"
    @info "Reduction blocks: $reduction_blocks"
    mod = create_module_specifics(model_param["Type"], module_list, @__MODULE__,
                                  "model_reduction_name")
    nmodes = get(model_param, "Number of Modes", 1)
    material_point_region = get(model_param, "Material Point Region", true)
    # Max Memory MiB: an OS/scheduler OOM kill is a SIGKILL, never catchable from inside
    # this process (see reduce_matrices' docstring); this lets the cascade stop itself
    # deliberately, with an ordinary Julia error naming the level, before that happens.
    # Defaults to 8192 (8 GiB), matching reduce_matrices' own default; an explicit
    # `Max Memory MiB: null` in the deck disables it.
    max_rss_mib = get(model_param, "Max Memory MiB", 8192.0)
    # Dense Shell Limit: above this many physical degrees of freedom, a cascade level's
    # own local system is built and factorized sparse instead of dense (see
    # reduce_matrices' docstring) -- matters only once a shell gets wide, so it has no
    # effect on the single-level schemes.
    dense_shell_limit = get(model_param, "Dense Shell Limit", 1500)
    extra_kwargs = model_param["Type"] == "Craig Bampton Cascade" ?
                   (;
                    max_rss_mib = isnothing(max_rss_mib) ? nothing :
                                  Float64(max_rss_mib),
                    dense_shell_limit = Int64(dense_shell_limit)) : (;)

    retained_nodes, condensed_nodes, pd_nodes,
    coupling_nodes = partition_nodes(block_nodes, reduction_blocks, material_point_region)
    @info "Model Reduction: $(length(retained_nodes)) retained, " *
          "$(length(condensed_nodes)) condensed, $(length(pd_nodes)) material point, " *
          "$(length(coupling_nodes)) coupling nodes"

    if pd_nodes != []
        # init_matrix already assembled every node, PD included, since it runs before
        # the reduction blocks (and therefore pd_nodes) are known. Rebuild from a blank
        # slate with them excluded now that they are known, rather than patch the
        # existing assembly: compute_model's subtract-old/add-new is a no-op here (old
        # and current state are still identical at this point in setup), so it would
        # leave their init_matrix contribution in place instead of removing it.
        nodes = setdiff(collect(1:Data_Manager.get_nnodes()), pd_nodes)
        @timeit "update_material_point_part" rebuild_matrix!(nodes)
        @info "Model Reduction: rebuilt stiffness matrix excluding material point nodes"
    end
    K = Data_Manager.get_stiffness_matrix()

    if retained_nodes == []
        @warn "No retained nodes defined for model reduction. Using full stiffness matrix."
        return
    end
    if condensed_nodes == []
        @warn "No condensed nodes defined for model reduction. Using full stiffness matrix."
        return
    end

    mark_coupling_field!(retained_nodes, condensed_nodes, pd_nodes, coupling_nodes)

    dof = Data_Manager.get_dof()
    nnodes = Data_Manager.get_nnodes()
    perm_retained = create_permutation(retained_nodes, dof, nnodes)
    perm_condensed = create_permutation(condensed_nodes, dof, nnodes)
    density_mass = expand_density_per_dof(density, dof)
    @info "Model Reduction: entering $(model_param["Type"]) with " *
          "$(length(perm_retained)) retained and $(length(perm_condensed)) condensed " *
          "degrees of freedom (dof=$dof)"

    # Reduced Matrix Cache: the cascade -- factorizing the condensed region's fixed-
    # interface eigenproblem, stage by stage -- is the part of setup this whole model
    # reduction feature exists to make affordable at all, but it still has to run once
    # per solver setup regardless. A parameter study that only ever changes the *free*
    # (non-reduced) region, never the reduced block's own material/geometry/reduction
    # settings, reruns the identical cascade for no reason; the paper this cascade
    # implements notes explicitly that its superelements "can be stored and reused...
    # when only the free region is changed". `retained_nodes`/`dof` are compared against
    # the cache on load as a minimal sanity check, not a full guarantee the deck is
    # otherwise unchanged -- an incompatible cache for the same path is the caller's own
    # responsibility to remove.
    cache_file = get(model_param, "Reduced Matrix Cache", nothing)
    loaded_from_cache = false
    if !isnothing(cache_file) && isfile(cache_file)
        @info "Model Reduction: loading cached reduced matrices from $cache_file"
        cached = Serialization.deserialize(cache_file)
        if cached.retained_nodes == retained_nodes && cached.dof == dof
            K_reduced = cached.K_reduced
            mass_reduced = cached.mass_reduced
            loaded_from_cache = true
        else
            @warn "Model Reduction: cached reduced matrices at $cache_file do not " *
                  "match this run's retained nodes/dof -- ignoring the cache and " *
                  "recomputing."
        end
    end
    if !loaded_from_cache
        @timeit "Condensation" K_reduced,
                               mass_reduced=mod.reduce_matrices(K, density_mass,
                                                                perm_retained,
                                                                perm_condensed, nmodes;
                                                                extra_kwargs...)
        if !isnothing(cache_file)
            @info "Model Reduction: saving reduced matrices to $cache_file"
            Serialization.serialize(cache_file,
                                    (K_reduced = K_reduced, mass_reduced = mass_reduced,
                                     retained_nodes = retained_nodes, dof = dof))
        end
    end

    dropzeros!(mass_reduced)
    dropzeros!(K_reduced)

    Data_Manager.set_stiffness_matrix(K_reduced)
    Data_Manager.set_mass_matrix(mass_reduced)

    Data_Manager.set_reduced_model_pd(pd_nodes)
    Data_Manager.set_reduced_model_retained(retained_nodes)

    @info "Model reduction is applied"
    @info "condensed $(length(condensed_nodes)), coupling $(length(coupling_nodes)), " *
          "material point $(length(pd_nodes))"
    return
end

end
