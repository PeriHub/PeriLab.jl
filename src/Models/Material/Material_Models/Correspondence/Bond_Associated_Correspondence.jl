# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Associated_Correspondence
#using LinearAlgebra
using StaticArrays
include("../../../../Support/helpers.jl")
include("../../../../Support/geometry.jl")
include("../../material_basis.jl")
using .Helpers: find_local_neighbors, invert, rotate, fastdot
using .Geometry:
    compute_strain,
    compute_bond_level_rotation_tensor,
    compute_bond_level_deformation_gradient
include("../../../Pre_calculation/pre_bond_associated_correspondence.jl")
using .Pre_Bond_Associated_Correspondence: compute_weighted_volume
using TimerOutputs

export init_model
export compute_model

"""
    correspondence_name()

Gives the correspondence material name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the material.

Example:
```julia
println(correspondence_name())
"Material Template"
```
"""
function correspondence_name()
    return "Bond Associated Correspondence"
end

function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    material_parameter::Dict,
)
    if !haskey(material_parameter, "Symmetry")
        @error "Symmetry for correspondence material is missing; options are 'isotropic plane strain', 'isotropic plane stress', 'anisotropic plane stress', 'anisotropic plane stress','isotropic' and 'anisotropic'. For 3D the plane stress or plane strain option is ignored."
        return nothing
    end
    if haskey(material_parameter, "Accuracy Order")
        datamanager.set_accuracy_order(material_parameter["Accuracy Order"])
    end

    dof = datamanager.get_dof()
    datamanager.create_bond_field("Bond Strain", Float64, "Matrix", dof)
    datamanager.create_bond_field("Bond Cauchy Stress", Float64, "Matrix", dof)
    datamanager.create_constant_bond_field("Bond Strain Increment", Float64, "Matrix", dof)
    datamanager.create_constant_node_field("Integral Nodal Stress", Float64, "Matrix", dof)

    datamanager.create_bond_field("Bond Rotation Tensor", Float64, "Matrix", dof)

    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    weighted_volume = datamanager.create_constant_node_field("Weighted Volume", Float64, 1)

    weighted_volume =
        compute_weighted_volume(nodes, nlist, volume, bond_damage, omega, weighted_volume)

    return datamanager
end

"""
    fields_to_local_synchronize()

Returns a user developer defined local synchronization. This happens before each model.

The structure of the Dict must because

    synchfield = Dict(
        "Field name" =>
            Dict("upload_to_cores" => true, "dof" => datamanager.get_dof()),
    )

or

    synchfield = Dict(
        "Field name" =>
            Dict("download_from_cores" => true, "dof" => datamanager.get_dof()),
    )

# Arguments

"""
function fields_to_local_synchronize()
    return Dict()
end

"""
https://link.springer.com/article/10.1007/s10409-021-01055-5


- Bond associated neighborhood is the overlap between nlist[iID] and nlist[nlist[iID][jID]]
- Filter equal nodes and create a new neighborhoodlist for bond -> bond_nlist
- calculate K, Kinv and defGrad -> already there if the neighborhood loop is in a function
- weighted volume (sum(volume(bond_nlist))/sum(volume[nlist[iID]]))
- global local IDs to be checked
  -> all neighbors search for neighbors at each core
  -> numbers are correct and it allows a change in size -> local ID is correct
"""

function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    material_parameter::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)

    rotation::Bool = datamanager.get_rotation()

    dof = datamanager.get_dof()
    nlist = datamanager.get_field("Neighborhoodlist")

    bond_damage = datamanager.get_bond_damage("NP1")
    horizon = datamanager.get_field("Horizon")
    omega = datamanager.get_field("Influence Function")
    volume = datamanager.get_field("Volume")

    bond_length = datamanager.get_field("Bond Length")

    bond_geometry = datamanager.get_field("Bond Geometry")
    bond_length = datamanager.get_field("Bond Length")
    bond_deformation = datamanager.get_field("Deformed Bond Geometry", "NP1")

    strain_N = datamanager.get_field("Bond Strain", "N")
    strain_NP1 = datamanager.get_field("Bond Strain", "NP1")

    stress_integral = datamanager.get_field("Integral Nodal Stress")
    cauchy_stress_N = datamanager.get_field("Cauchy Stress", "N")
    cauchy_stress_NP1 = datamanager.get_field("Cauchy Stress", "NP1")
    stress_N = datamanager.get_field("Bond Cauchy Stress", "N")
    stress_NP1 = datamanager.get_field("Bond Cauchy Stress", "NP1")
    strain_increment_nodal = datamanager.get_field("Strain Increment")
    strain_increment = datamanager.get_field("Bond Strain Increment")
    bond_force = datamanager.get_field("Bond Forces")

    # computed in pre calculation ----------------------------------
    gradient_weights = datamanager.get_field("Lagrangian Gradient Weights")
    weighted_volume = datamanager.get_field("Weighted Volume")
    deformation_gradient = datamanager.get_field("Weighted Deformation Gradient")
    #---------------------------------------------------------------
    displacements = datamanager.get_field("Displacements", "NP1")
    velocity = datamanager.get_field("Velocity", "NP1")

    ba_deformation_gradient = datamanager.get_field("Bond Associated Deformation Gradient")

    ba_deformation_gradient = compute_bond_level_deformation_gradient(
        nodes,
        nlist,
        dof,
        bond_geometry,
        bond_length,
        bond_deformation,
        deformation_gradient,
        ba_deformation_gradient,
    )

    ba_rotation_tensor = datamanager.get_field("Bond Rotation Tensor", "NP1")

    strain_NP1 = compute_bond_strain(nodes, nlist, ba_deformation_gradient, strain_NP1)
    strain_increment = strain_NP1 - strain_N
    # TODO decomposition to get the rotation and large deformation in
    # TODO store not angles, but rotation matrices, because they are computed in decomposition
    if rotation
        rotation_tensor = datamanager.get_field("Rotation Tensor")
        ba_rotation_tensor = compute_bond_level_rotation_tensor(
            nodes,
            nlist,
            ba_deformation_gradient,
            ba_rotation_tensor,
        )
        nneighbors = datamanager.get_field("Number of Neighbors")
        for iID in nodes
            stress_N[iID] = rotate(
                Vector{Int64}(1:nneighbors[iID]),
                stress_N[iID],
                ba_rotation_tensor[iID],
                false,
            )
            strain_increment[iID] = rotate(
                Vector{Int64}(1:nneighbors[iID]),
                strain_increment[iID],
                ba_rotation_tensor[iID],
                false,
            )
        end
    end

    material_models = split(material_parameter["Material Model"], "+")
    material_models = map(r -> strip(r), material_models)

    for material_model in material_models
        mod = datamanager.get_model_module(material_model)
        for iID in nodes

            #cauchy_stress_NP1, datamanager = mod.compute_stresses(datamanager, iID, dof, material_parameter, time, dt, strain_increment_nodal, cauchy_stress_N, stress_NP1)

            for (jID, nID) in enumerate(nlist[iID])
                # TODO how to make the separation if the datamager is included?
                stress_NP1[iID], datamanager = mod.compute_stresses(
                    datamanager,
                    jID,
                    dof,
                    material_parameter,
                    time,
                    dt,
                    strain_increment[iID],
                    stress_N[iID],
                    stress_NP1[iID],
                    (iID, jID, nID),
                )
            end
        end
    end
    if rotation
        for iID in nodes
            stress_NP1[iID] = rotate(
                Vector{Int64}(1:nneighbors[iID]),
                stress_NP1[iID],
                ba_rotation_tensor[iID],
                true,
            )
        end
    end

    stress_integral = compute_stress_integral(
        nodes,
        dof,
        nlist,
        omega,
        bond_damage,
        volume,
        weighted_volume,
        bond_geometry,
        bond_length,
        stress_NP1,
        ba_deformation_gradient,
        stress_integral,
    )

    bond_force = compute_bond_forces(
        nodes,
        nlist,
        bond_geometry,
        bond_length,
        stress_NP1,
        stress_integral,
        weighted_volume,
        gradient_weights,
        omega,
        bond_damage,
        bond_force,
    )

    return datamanager

end

function compute_stress_integral(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    nlist::Union{Vector{Vector{Int64}},SubArray},
    omega::SubArray,
    bond_damage::SubArray,
    volume::SubArray,
    weighted_volume::SubArray,
    bond_geometry::SubArray,
    bond_length::SubArray,
    bond_stresses::SubArray,
    deformation_gradient::SubArray,
    stress_integral::SubArray,
)
    temp = @MMatrix zeros(dof, dof)
    one = @views I(dof)
    factor::Float64 = 0
    for iID in nodes
        stress_integral[iID, :, :] .= 0.0
        for (jID, nID) in enumerate(nlist[iID])
            if bond_damage[iID][jID] == 0
                continue
            end
            temp =
                (one - bond_geometry[iID][jID, :] * bond_geometry[iID][jID, :]') ./
                (bond_length[iID][jID] * bond_length[iID][jID])
            factor =
                volume[nID] *
                omega[iID][jID] *
                bond_damage[iID][jID] *
                (0.5 / weighted_volume[iID] + 0.5 / weighted_volume[nID])

            stress_integral[iID, :, :] +=
                factor .* compute_Piola_Kirchhoff_stress(
                    bond_stresses[iID][jID, :, :],
                    deformation_gradient[iID][jID, :, :],
                ) * temp

        end
    end
    return stress_integral
end

#function compute_bond_strain(nodes::Union{SubArray,Vector{Int64}}, nlist::Union{Vector{Vector{Int64}},SubArray}, deformation_gradient::SubArray, strain::SubArray)
#
function compute_bond_strain(nodes, nlist, deformation_gradient, strain)

    for iID in nodes
        strain[iID][:, :, :] = compute_strain(
            eachindex(nlist[iID]),
            (@view deformation_gradient[iID][:, :, :]),
            (@view strain[iID][:, :, :]),
        )
    end
    return strain
end





function update_Green_Langrange_nodal_strain_increment(
    nodes::Union{SubArray,Vector{Int64}},
    dt::Float64,
    deformation_gradient::SubArray,
    deformation_gradient_dot::SubArray,
    strain_increment::SubArray,
)

    for iID in nodes
        strain_increment[iID, :, :] = update_Green_Langrange_strain(
            dt,
            deformation_gradient[iID, :, :],
            deformation_gradient_dot[iID, :, :],
            strain_increment[iID, :, :],
        )
    end

end

function update_Green_Langrange_bond_strain_increment(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::Union{Vector{Vector{Int64}},SubArray},
    dt::Float64,
    deformation_gradient::SubArray,
    deformation_gradient_dot::SubArray,
    strain_increment::SubArray,
)

    for iID in nodes
        for jID in nlist[iID]
            update_Green_Langrange_strain(
                dt,
                deformation_gradient[iID][jID, :, :],
                deformation_gradient_dot[iID][jID, :, :],
                strain_increment[iID][jID, :, :],
            )
        end
    end

end


function update_Green_Langrange_strain(
    dt::Float64,
    deformation_gradient::Matrix{Float64},
    deformation_gradient_dot::Matrix{Float64},
    strain::Matrix{Float64},
)

    @inbounds @fastmath for m ∈ axes(deformation_gradient, 1),
        n ∈ axes(deformation_gradient_dot, 2)

        strain_mn = zero(Float64)
        for k ∈ axes(deformation_gradient, 2)
            strain_mn += deformation_gradient[m, k] * deformation_gradient_dot[k, n]
        end
        strain[m, n] = strain_mn * dt * 0.5
        strain[n, m] = strain_mn
    end

end

function compute_bond_forces(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::Union{Vector{Vector{Int64}},SubArray},
    bond_geometry,
    bond_length,
    bond_stress,
    integral_nodal_stress,
    weighted_volume,
    gradient_weights,
    omega,
    bond_damage,
    bond_forces,
)

    for iID in nodes

        for (jID, nID) in enumerate(nlist[iID])
            if bond_damage[iID][jID] == 0
                continue
            end


            bond_forces[iID][jID, :] =
                integral_nodal_stress[iID, :, :] * gradient_weights[iID][jID, :]

            bond_forces[iID][jID, :] +=
                bond_damage[iID][jID] * omega[iID][jID] /
                (weighted_volume[iID] * bond_length[iID][jID] * bond_length[iID][jID]) .*
                bond_stress[iID][jID, :, :] * bond_geometry[iID][jID, :]


            #
            #mul!(bond_forces[iID][jID, :], integral_nodal_stress[iID, :, :], gradient_weights[iID][jID, :])
        end
    end
    return bond_forces
end


end
