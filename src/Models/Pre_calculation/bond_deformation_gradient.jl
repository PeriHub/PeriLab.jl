# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation_Gradient

include("../../Support/geometry.jl")
using .Geometry: calculate_bond_length, compute_weighted_deformation_gradient
export pre_calculation_name
export init_model
export compute

include("../../Support/helpers.jl")
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
"Bond Deformation Gradient"
```
"""
function pre_calculation_name()
    return "Bond Deformation Gradient"
end


"""
    init_model(datamanager, nodes, parameter)

Inits the bond-based corrosion model. This template has to be copied, the file renamed and edited by the user to create a new corrosion. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    parameter::Dict,
)
    dof = datamanager.get_dof()
    datamanager.create_constant_bond_field(
        "Bond Associated Deformation Gradient",
        Float64,
        "Matrix",
        dof,
    )
    datamanager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    datamanager.create_constant_bond_field("Lagrangian Gradient Weights", Float64, dof)
    datamanager.create_constant_node_field(
        "Weighted Deformation Gradient",
        Float64,
        "Matrix",
        dof,
    )
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


function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

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

    weighted_volume =
        compute_weighted_volume(nodes, nlist, volume, bond_damage, omega, weighted_volume)
    gradient_weights = compute_Lagrangian_gradient_weights(
        nodes,
        dof,
        accuracy_order,
        volume,
        nlist,
        horizon,
        bond_damage,
        omega,
        bond_geometry,
        gradient_weights,
    )
    deformation_gradient = compute_weighted_deformation_gradient(
        nodes,
        dof,
        nlist,
        volume,
        gradient_weights,
        displacement,
        deformation_gradient,
    )
    return datamanager
end



function compute_weighted_volume(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::Union{Vector{Vector{Int64}},SubArray},
    volume::SubArray,
    bond_damage::SubArray,
    omega::SubArray,
    weighted_volume::SubArray,
)
    for iID in nodes
        weighted_volume[iID] =
            sum(bond_damage[iID][:] .* omega[iID][:] .* volume[nlist[iID]])
    end
    return weighted_volume
end


"""
accuracy_order::Int64 - needs a number of bonds which are linear independent

"""

function calculate_Q(
    accuracy_order::Int64,
    dof::Int64,
    bond_geometry::Vector{Float64},
    horizon::Union{Int64,Float64},
)

    Q = ones(Float64, qdim(accuracy_order, dof))  # Initialize Q with ones
    counter = 0
    p = zeros(Int64, dof)
    for this_order = 1:accuracy_order
        for p[1] = this_order:-1:0
            if dof == 3
                for p[2] = this_order-p[1]:-1:0
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
function compute_Lagrangian_gradient_weights(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    accuracy_order::Int64,
    volume::Union{SubArray,Vector{Float64}},
    nlist::Union{Vector{Vector{Int64}},SubArray},
    horizon::Union{SubArray,Vector{Float64}},
    bond_damage::Union{SubArray,Vector{Vector{Float64}}},
    omega::Union{SubArray,Vector{Vector{Float64}}},
    bond_geometry,
    gradient_weights,
)
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

        Minv = invert(
            M,
            "In compute_Lagrangian_gradient_weights the matrix M is singular and cannot be inverted. To many bond damages or a to small horizon might cause this.",
        )

        for (jID, nID) in enumerate(nlist[iID])
            Q = calculate_Q(accuracy_order, dof, bond_geometry[iID][jID, :], horizon[iID])
            # this comes from Eq(19) in 10.1007/s40571-019-00266-9
            # or example 1 in https://arxiv.org/pdf/2004.11477
            for idof = 1:dof # Eq (3) flowing
                gradient_weights[iID][jID, idof] =
                    omega[iID][jID] * bond_damage[iID][jID] / horizon[iID] *
                    Minv[idof, :]' *
                    Q
            end

        end
    end
    return gradient_weights
end

end
