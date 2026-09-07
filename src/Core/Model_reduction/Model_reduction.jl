# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Model_reduction
using TimerOutputs: @timeit
using ...Data_Manager
using SparseArrays
using ...Solver_Manager: find_module_files, create_module_specifics
using ...Correspondence_matrix_based: build_mass_matrix, init_model, init_matrix,
                                      compute_model
using ...Helpers: find_active_nodes, create_permutation
global module_list = find_module_files(@__DIR__, "model_reduction_name")
for mod in module_list
    include(mod["File"])
end

export ReducedState
export n_physical, n_modal, physical_part, modal_part
export pull_from_nodes!, push_to_nodes!
export setup_reduced_state

"""
    ReducedState

State vectors of the time integration, held flat rather than as node fields.

The vector is ordered `[u_master; eta]`: `n_phys` physical degrees of freedom of the
master nodes, followed by `n_modal` modal coordinates. A modal coordinate belongs to no
node and therefore cannot live in a `nnodes x dof` field, which is what this container
is for.

`n_modal` is zero for a static reduction (Guyan) and for a run without any reduction, in
which case the state is nothing but the master degrees of freedom and the integration is
the same as before. The unreduced case is the one where every node is a master.

The physical part is ordered component by component, the same way `vec` flattens an
`n x dof` node field: first the first component of every master node, then the second,
and so on. Master node `k` therefore sits at `(d - 1) * n_master + k` for component `d`.
This has to match the ordering the reduced matrices were built in.

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

View on the entries that belong to master nodes.

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
    pull_from_nodes!(vector, field, master_nodes, state)

Copies the master rows of a node field into the physical part of a state vector.

The modal part is left untouched, because a node field says nothing about it.

# Arguments
- `vector::AbstractVector{Float64}`: Target, length `n_phys + n_modal`
- `field::AbstractMatrix{Float64}`: Node field, `nnodes x dof`
- `master_nodes::AbstractVector{Int64}`: Master node indices
- `state::ReducedState`: The state, for the split position
# Returns
- `vector`: The unchanged reference
"""
function pull_from_nodes!(vector::AbstractVector{Float64},
                          field::AbstractMatrix{Float64},
                          master_nodes::AbstractVector{Int64},
                          state::ReducedState)
    dof = size(field, 2)
    n_master = length(master_nodes)
    # Component major, matching vec(field[master_nodes, :]).
    @inbounds for d in 1:dof
        offset = (d - 1) * n_master
        for (k, node) in enumerate(master_nodes)
            vector[offset + k] = field[node, d]
        end
    end
    return vector
end

"""
    push_to_nodes!(field, vector, master_nodes, state)

Copies the physical part of a state vector back into the master rows of a node field.

Rows of nodes that are no master are left alone, and the modal part is dropped: it has
no node to be written to.

# Arguments
- `field::AbstractMatrix{Float64}`: Node field, `nnodes x dof`
- `vector::AbstractVector{Float64}`: Source, length `n_phys + n_modal`
- `master_nodes::AbstractVector{Int64}`: Master node indices
- `state::ReducedState`: The state, for the split position
# Returns
- `field`: The unchanged reference
"""
function push_to_nodes!(field::AbstractMatrix{Float64},
                        vector::AbstractVector{Float64},
                        master_nodes::AbstractVector{Int64},
                        state::ReducedState)
    dof = size(field, 2)
    n_master = length(master_nodes)
    # Component major, matching vec(field[master_nodes, :]).
    @inbounds for d in 1:dof
        offset = (d - 1) * n_master
        for (k, node) in enumerate(master_nodes)
            field[node, d] = vector[offset + k]
        end
    end
    return field
end

"""
    setup_reduced_state(model_reduction, K)

Master nodes and state vectors for the time loop.

Without a reduction every node is a master and there are no modal coordinates, so the
state is the full system written as a flat vector and the time loop needs no branch on
whether a reduction is active.

# Arguments
- `model_reduction`: The `"Model Reduction"` solver option, `false` when disabled
- `K::AbstractMatrix`: The stiffness matrix, reduced or not
# Returns
- `master_nodes::Vector{Int64}`: Nodes carrying physical degrees of freedom
- `state::ReducedState`: Zero initialised state
"""
function setup_reduced_state(model_reduction, K::AbstractMatrix)
    dof = Data_Manager.get_dof()

    master_nodes = model_reduction == false ?
                   collect(1:Data_Manager.get_nnodes()) :
                   Data_Manager.get_reduced_model_master()

    n_phys = length(master_nodes) * dof
    n_total = size(K, 1)

    if n_total < n_phys
        throw(ArgumentError("Stiffness matrix has $n_total degrees of freedom for " *
                            "$n_phys master degrees of freedom. Master node list and " *
                            "matrix disagree."))
    end

    state = ReducedState(n_phys, n_total)
    if n_modal(state) > 0
        @info "Reduced system: $n_phys physical and $(n_modal(state)) modal degrees of freedom"
    end

    return master_nodes, state
end

function reduce_model(K::AbstractMatrix{Float64}, M::AbstractMatrix{Float64},
                      m::Vector{Int64}, s::Vector{Int64}; n_modes::Int64 = 1)
    mod = Data_Manager.get_model_module(model_param["Reduction Model"])

    mod.reduced_matrices(K, M, m, s; n_modes)

    if isnothing(Data_Manager.get_filtered_nlist())
        @timeit "compute index" return damage_index(nodes)
    end
end

function init_reduce_model(model_param::Dict, block_nodes::Dict{Int64,Vector{Int64}},
                           density)
    reduction_blocks = []

    reduction_blocks = get(model_param, "Reduction Blocks", nothing)

    if isnothing(reduction_blocks)
        @warn "No reduction blocks defined for model reduction. If you want to use a reduced model please define 'Reduction Blocks' in the yaml input deck."
        return
    end
    if reduction_blocks isa Float64
        @error "Type Float is not supported for Reduction Blocks"
        return
    end

    if reduction_blocks isa Int64
        reduction_blocks = [reduction_blocks]
    else
        reduction_blocks = parse.(Int64,
                                  filter(!isempty,
                                         split(model_param["Reduction Blocks"],
                                               r"[,\s]")))
    end

    @info "Model Reduction Type: $(model_param["Type"])"
    mod = create_module_specifics(model_param["Type"],
                                  module_list,
                                  @__MODULE__,
                                  "model_reduction_name")

    nmodes = get(model_param, "Number of Modes", 1)
    master_nodes = Int64[]
    slave_nodes = Int64[]
    pd_nodes = Int64[]

    nlist = Data_Manager.get_nlist()
    # only for visualization and debugging.
    cn = Data_Manager.create_constant_node_scalar_field("Coupling Nodes", Int64)
    full_blocks = setdiff(collect(keys(block_nodes)), reduction_blocks)
    @info "Reduction blocks: $reduction_blocks"
    for block in full_blocks
        append!(pd_nodes, block_nodes[block])
    end

    for node in pd_nodes
        append!(master_nodes, nlist[node])
    end

    append!(master_nodes, pd_nodes)

    master_nodes = sort(unique(master_nodes))

    if !(get(model_param, "Material Point Region", true))
        pd_nodes::Vector{Int64} = []
    end
    sort!(pd_nodes)
    slave_nodes = sort!(setdiff(collect(1:Data_Manager.get_nnodes()), master_nodes))
    if pd_nodes != []
        nodes = setdiff(collect(1:Data_Manager.get_nnodes()), pd_nodes)
        # update matrix excluding PD nodes
        @timeit "update_material_point_part" compute_model(nodes)
    end
    K = Data_Manager.get_stiffness_matrix()
    if master_nodes == []
        @warn "No master nodes defined for model reduction. Using full stiffness matrix."
        return
    end
    if slave_nodes == []
        @warn "No slave nodes defined for model reduction. Using full stiffness matrix."
        return
    end
    coupling_nodes = setdiff(master_nodes, pd_nodes)

    cn[master_nodes] .= 3
    cn[slave_nodes] .= 6
    cn[pd_nodes] .+= 1  # added to be sure that all points are handled
    cn[coupling_nodes] .+= 3  # added to be sure that all points are handled

    nnodes = Data_Manager.get_nnodes()
    dof = Data_Manager.get_dof()
    perm_master = create_permutation(master_nodes, Data_Manager.get_dof(), nnodes)
    perm_slave = create_permutation(slave_nodes, Data_Manager.get_dof(), nnodes)

    # create the mass part.

    density_mass = zeros(Float64, length(density) * dof)

    for iID in eachindex(density_mass)
        density_mass[iID] = density[Int(ceil(iID / dof))]
    end
    # perform the condensation of the system

    @timeit "Condensation" K_reduced,
                           mass_reduced=mod.reduce_matrices(K,
                                                            density_mass,
                                                            perm_master,
                                                            perm_slave,
                                                            nmodes)

    dropzeros!(mass_reduced)
    dropzeros!(K_reduced)

    Data_Manager.set_stiffness_matrix(K_reduced)
    Data_Manager.set_mass_matrix(mass_reduced)

    Data_Manager.set_reduced_model_pd(pd_nodes)
    Data_Manager.set_reduced_model_master(master_nodes)

    @info "Model reduction is applied"
    @info "Condensed: $(length(slave_nodes)), Master: $(length(master_nodes)), PD: $(length(pd_nodes))."
    return
end

end
