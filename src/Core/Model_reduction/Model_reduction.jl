# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Model_reduction
using TimerOutputs: @timeit
using ....Data_Manager
using ...Solver_Manager: find_module_files, create_module_specifics
using ...Correspondence_matrix_based: build_mass_matrix, init_model, init_matrix,
                                      compute_model
using ...Helpers: find_active_nodes, create_permutation
global module_list = find_module_files(@__DIR__, "model_reduction_name")
for mod in module_list
    include(mod["File"])
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

    if pd_nodes != []
        nodes = setdiff(collect(1:Data_Manager.get_nnodes()), pd_nodes)
        # update matrix excluding PD nodes
        @timeit "update_material_point_part" compute_model(nodes)
    end
    K = Data_Manager.get_stiffness_matrix()
    if !(master_nodes == [])
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

        density_mass = zeros(length(density) * dof)

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
        @info "b"

        dropzeros!(mass_reduced)
        dropzeros!(K_reduced)

        Data_Manager.set_stiffness_matrix(sparse(K_reduced))
        Data_Manager.set_mass_matrix(sparse(mass_reduced))

        Data_Manager.set_reduced_model_pd(pd_nodes)
        Data_Manager.set_reduced_model_master(master_nodes)

        @info "Model reduction is applied"
        @info "Condensed: $(length(slave_nodes)), Master: $(length(master_nodes)), PD: $(length(pd_nodes))."
        return
    end
    @warn "No master nodes defined for model reduction. Using full stiffness matrix."
end

end
