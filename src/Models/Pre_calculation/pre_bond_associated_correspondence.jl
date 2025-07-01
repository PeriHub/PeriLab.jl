# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_Bond_Associated_Correspondence
using DataStructures: OrderedDict
include("../../Support/Geometry.jl")
using .Geometry: compute_weighted_deformation_gradient
using LoopVectorization
using StaticArrays: @MVector
export fields_for_local_synchronization
export pre_calculation_name
export init_model
export compute

include("../../Support/Helpers.jl")
using .Helpers: invert, qdim

"""
    pre_calculation_name()

Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the Pre_Calculation.

Example:
```julia
println(pre_calculation_name())
"Bond Associated Correspondence"
```
"""
function pre_calculation_name()
    return "Bond Associated Correspondence"
end

"""
    init_model(datamanager, nodes, parameter)

Inits the bond deformation gradient calculation.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    parameter::Union{Dict,OrderedDict},
                    block::Int64)
    dof = datamanager.get_dof()
    datamanager.create_constant_bond_field("Bond Associated Deformation Gradient",
                                           Float64,
                                           "Matrix",
                                           dof)
    datamanager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    datamanager.create_constant_bond_field("Lagrangian Gradient Weights", Float64, dof)
    datamanager.create_constant_node_field("Weighted Deformation Gradient",
                                           Float64,
                                           "Matrix",
                                           dof)

    #https://arxiv.org/pdf/2004.11477
    # maybe as static array

    accuracy_order = datamanager.get_accuracy_order()
    dim = qdim(accuracy_order, dof)
    if "Q Vector" in datamanager.get_all_field_keys()
        return datamanager
    end
    datamanager.create_constant_free_size_field("Q Vector", Float64, (dim,), 1)
    datamanager.create_constant_free_size_field("Minv Matrix", Float64, (dim, dim))
    datamanager.create_constant_free_size_field("M temporary Matrix", Float64, (dim, dim))
    return datamanager
end

"""
    compute(datamanager, nodes)

Compute the bond deformation gradient.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
# Returns
- `datamanager`: Datamanager.
"""

function compute(datamanager::Module,
                 nodes::AbstractVector{Int64},
                 parameter::Union{Dict,OrderedDict},
                 block::Int64)
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    bond_geometry = datamanager.get_field("Bond Geometry")
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    displacement = datamanager.get_field("Displacements", "NP1")
    gradient_weights = datamanager.get_field("Lagrangian Gradient Weights")
    weighted_volume = datamanager.get_field("Weighted Volume")
    bond_length = datamanager.get_field("Bond Length")
    deformation_gradient = datamanager.get_field("Weighted Deformation Gradient")
    horizon = datamanager.get_field("Horizon")
    ba_deformation_gradient = datamanager.get_field("Bond Associated Deformation Gradient")
    ba_rotation_tensor = datamanager.get_field("Bond Rotation Tensor", "NP1")
    accuracy_order = datamanager.get_accuracy_order()

    Q = datamanager.get_field("Q Vector")
    Minv = datamanager.get_field("Minv Matrix")
    M = datamanager.get_field("M temporary Matrix")

    compute_weighted_volume!(weighted_volume, nodes, nlist, volume, bond_damage, omega)

    compute_Lagrangian_gradient_weights(nodes,
                                        dof,
                                        accuracy_order,
                                        volume,
                                        nlist,
                                        horizon,
                                        bond_damage,
                                        omega,
                                        Q,
                                        M,
                                        Minv,
                                        bond_geometry,
                                        gradient_weights)
    compute_weighted_deformation_gradient(nodes,
                                          nlist,
                                          volume,
                                          gradient_weights,
                                          displacement,
                                          deformation_gradient)

    return datamanager
end

function compute_weighted_volume!(weighted_volume::Vector{Float64},
                                  nodes::AbstractVector{Int64},
                                  nlist::Union{Vector{Vector{Int64}},SubArray},
                                  volume::Vector{Float64},
                                  bond_damage::Vector{Vector{Float64}},
                                  omega::Vector{Vector{Float64}})
    for iID in nodes
        weighted_volume[iID] = 0
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            weighted_volume[iID] += bond_damage[iID][jID] * omega[iID][jID] *
                                    volume[nlist[iID][jID]]
        end
    end
end

"""
accuracy_order::Int64 - needs a number of bonds which are linear independent

"""

function calculate_Q_old(accuracy_order::Int64,
                         dof::Int64,
                         bond_geometry::Vector{Float64},
                         horizon::Union{Int64,Float64})
    Q = ones(Float64, qdim(accuracy_order, dof))  # Initialize Q with ones
    counter = 0
    p = zeros(Int64, dof)
    for this_order in 1:accuracy_order
        for p[1] in this_order:-1:0
            if dof == 3
                for p[2] in (this_order - p[1]):-1:0
                    p[3] = this_order - p[1] - p[2]
                    # Calculate the product for Q[counter]
                    counter += 1
                    Q[counter] = prod((bond_geometry ./ horizon) .^ p)
                end
            else
                p[2] = this_order - p[1]
                counter += 1
                Q[counter] = prod((bond_geometry ./ horizon) .^ p)
            end
        end
    end
    return Q
end

function calculate_Q(accuracy_order::Int64,
                     dof::Int64,
                     bond_geometry::Vector{Float64},
                     horizon::Union{Int64,Float64},
                     Q::Vector{Float64})
    counter = 0
    p = @MVector zeros(Int64, dof)
    for this_order in 1:accuracy_order
        for p[1] in this_order:-1:0
            if dof == 3
                for p[2] in (this_order - p[1]):-1:0
                    p[3] = this_order - p[1] - p[2]
                    # Calculate the product for Q[counter]
                    counter += 1
                    Q[counter] = prod_Q(bond_geometry, horizon, p, Q[counter])
                end
            else
                p[2] = this_order - p[1]
                counter += 1
                Q[counter] = prod_Q(bond_geometry, horizon, p, Q[counter])
            end
        end
    end
    return Q
end

function prod_Q(bond_geometry, horizon, p, Q)
    Q = 1
    @inbounds @fastmath for m in axes(p, 1)
        @views @inbounds @fastmath for pc in 1:p[m]
            Q *= bond_geometry[m] / horizon
        end
    end
    return Q
end

function compute_Lagrangian_gradient_weights(nodes::AbstractVector{Int64},
                                             dof::Int64,
                                             accuracy_order::Int64,
                                             volume::AbstractVector{Float64},
                                             nlist::Union{Vector{Vector{Int64}},SubArray},
                                             horizon::Union{SubArray,Vector{Float64}},
                                             bond_damage::Union{SubArray,
                                                                Vector{Vector{Float64}}},
                                             omega::Union{SubArray,
                                                          Vector{Vector{Float64}}},
                                             Q,
                                             M,
                                             Minv,
                                             bond_geometry,
                                             gradient_weights)
    for iID in nodes
        M .= 0
        # bg = @view bond_geometry[iID]
        for (jID, nID) in enumerate(nlist[iID])
            Q = calculate_Q(accuracy_order, dof, bond_geometry[iID][jID], horizon[iID], Q)
            QTQ!(M, omega[iID][jID], bond_damage[iID][jID], volume[nID], Q)
        end

        @views Minv = invert(M,
                             "In compute_Lagrangian_gradient_weights the matrix M is singular and cannot be inverted. To many bond damages or a to small horizon might cause this.")

        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            Q = calculate_Q(accuracy_order, dof, bond_geometry[iID][jID], horizon[iID], Q)
            # this comes from Eq(19) in 10.1007/s40571-019-00266-9
            # or example 1 in https://arxiv.org/pdf/2004.11477
            # Eq (3) flowing
            compute_gradient_weights!(gradient_weights,
                                      dof,
                                      omega,
                                      bond_damage,
                                      horizon,
                                      Minv,
                                      Q,
                                      iID,
                                      jID)
        end
    end
end
function QTQ!(M, omega, bond_damage, volume, Q)
    # each operation does not need the sum. It's an integral over all neighbor nodes.
    @inbounds @fastmath for m in axes(M, 1)
        @inbounds @fastmath for n in axes(M, 2)
            M[m, n] += omega * bond_damage * volume * Q[m] * Q[n]
        end
    end
end

function compute_gradient_weights!(gradient_weights,
                                   dof,
                                   omega,
                                   bond_damage,
                                   horizon,
                                   Minv,
                                   Q,
                                   iID,
                                   jID)
    for idof in 1:dof # Eq (3) flowing
        gradient_weights[iID][jID][idof] = 0
        @inbounds @fastmath for m in axes(Minv, 2)
            gradient_weights[iID][jID][idof] += omega[iID][jID] * bond_damage[iID][jID] /
                                                horizon[iID] *
                                                Minv[idof, m] *
                                                Q[m]
        end
    end
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    # download_from_cores = false
    # upload_to_cores = true
    # datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

end
