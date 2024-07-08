# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation_Gradient

include("../Material/Material_Models/Bond_Associated_Correspondence.jl")
using .Bond_Associated_Correspondence: find_local_neighbors
include("../../Support/geometry.jl")
using .Geometry: calculate_bond_length, compute_bond_level_deformation_gradient
export compute
include("../../Support/helpers.jl")
using .Helpers: invert, qdim

"""
    compute(datamanager, nodes)

Compute the bond deformation gradient.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
# Returns
- `datamanager`: Datamanager.
"""


function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block_id::Int64)

    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    gradient_weights = datamanager.get_field("Lagrangian Gradient Weights")

    gradient_weights = datamanager.get_field("Lagrangian Gradient Weights")

    weighted_volume = datamanager.get_field("Weighted Volume")

    displacements = datamanager.get_field("Displacements", "NP1")
    velocity = datamanager.get_field("Velocity", "NP1")
    accuracy_order = material_parameter["Accuracy Order"]

    deformation_gradient = datamanager.get_field("Deformation Gradient")

    ba_deformation_gradient = datamanager.get_field("Bond Deformation Gradient", "NP1")
    ba_rotation_tensor = datamanager.get_field("Bond Rotation Tensor", "NP1")
    try
        accuracy_order = datamanager.get_property(block_id, "Material Model", "Accuracy Order")
    catch
        @error "Accuracy order must be defined in all Materials if one is bond associated"
    end
    weighted_volume = compute_weighted_volume(nodes, nlist, volume, bond_damage, omega, weighted_volume)
    gradient_weights = compute_Lagrangian_gradient_weights(nodes, dof, accuracy_order, volume, nlist, horizon, bond_damage, omega, bond_geometry, gradient_weights)

    ba_deformation_gradient = compute_bond_level_deformation_gradient(nodes, nlist, dof, bond_geometry, bond_length, bond_deformation, deformation_gradient, ba_deformation_gradient)


    return datamanager
end



function compute_weighted_volume(nodes::Union{SubArray,Vector{Int64}}, nlist::Union{Vector{Vector{Int64}},SubArray}, volume::SubArray, bond_damage::SubArray, omega::SubArray, weighted_volume::SubArray)
    for iID in nodes
        weighted_volume[iID] = sum(bond_damage[iID][:] .* omega[iID][:] .* volume[nlist[iID]])
    end
    return weighted_volume
end


"""
accuracy_order::Int64 - needs a number of bonds which are linear independent

"""

function calculate_Q(accuracy_order::Int64, dof::Int64, bond_geometry::Vector{Float64}, horizon::Union{Int64,Float64})

    Q = ones(Float64, qdim(accuracy_order, dof))  # Initialize Q with ones
    counter = 1
    p = zeros(Int64, dof)
    for this_order in 1:accuracy_order
        for p[1] in this_order:-1:0
            if dof == 3
                for p[2] in this_order-p[1]:-1:0
                    p[3] = this_order - p[1] - p[2]
                    # Calculate the product for Q[counter]
                    Q[counter] = prod((bond_geometry ./ horizon) .^ p)
                    counter += 1
                end
            else
                p[2] = this_order - p[1]
                Q[counter] = prod((bond_geometry ./ horizon) .^ p)
                counter += 1
            end
        end
    end
    return Q
end
function compute_Lagrangian_gradient_weights(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, accuracy_order::Int64, volume::Union{SubArray,Vector{Float64}}, nlist::Union{Vector{Vector{Int64}},SubArray}, horizon::Union{SubArray,Vector{Float64}}, bond_damage::Union{SubArray,Vector{Vector{Float64}}}, omega::Union{SubArray,Vector{Vector{Float64}}}, bond_geometry, gradient_weights)
    #https://arxiv.org/pdf/2004.11477
    # maybe as static array
    dim = qdim(accuracy_order, dof)
    Minv = zeros(Float64, dim, dim)
    for iID in nodes
        M = zeros(Float64, dim, dim)
        for (jID, nID) in enumerate(nlist[iID])
            Q = calculate_Q(accuracy_order, dof, bond_geometry[iID][jID, :], horizon[iID])
            M += omega[iID][jID] * bond_damage[iID][jID] * volume[nID] .* Q * Q'
        end
        try
            Minv = inv(M)
        catch
            @error "In compute_Lagrangian_gradient_weights the matrix M is singular and cannot be inverted. To many bond damages or a to small horizon might cause this."
            return nothing
        end
        for (jID, nID) in enumerate(nlist[iID])
            Q = calculate_Q(accuracy_order, dof, bond_geometry[iID][jID, :], horizon[iID])
            # this comes from Eq(19) in 10.1007/s40571-019-00266-9
            # or example 1 in https://arxiv.org/pdf/2004.11477
            for idof in 1:dof
                gradient_weights[iID][jID, idof] = omega[iID][jID] / horizon[iID] .* (Minv[idof, :]' * Q)
            end
        end
    end
    return gradient_weights
end

end