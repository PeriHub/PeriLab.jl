# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Elastic
using .......Data_Manager
using .......ParameterSpec: @params, register_material
using .....Material_Basis: hooke_matrix
using .......Helpers: get_fourth_order, fast_mul!, get_mapping
using StaticArrays: SMatrix
export compute_stresses
export correspondence_name
export fe_support
export init_model
export fields_for_local_synchronization

"Parameters of Correspondence Elastic beyond the shared material keys (none)."
@params struct CorrespondenceElasticParams
end
__init__() = register_material("Correspondence Elastic", CorrespondenceElasticParams)

"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
    return true
end

"""
	correspondence_name()

Gives the material name. PeriLab loads the module because it defines this function; the input deck uses the name passed to `register_*` in `__init__()`.

# Arguments

# Returns
- `name::String`: The name of the material.

Example:
```julia
println(material_name())
"Material Template"
```
"""
function correspondence_name()
    return "Correspondence Elastic"
end

"""
	compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, p, material, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

Calculates the stresses of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `dof::Int64`: Degrees of freedom
- `p`: The model parameters.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{NodeTensorField{Float64},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
- `iID_jID_nID::Tuple=(): (optional) are the index and node id information. The tuple is ordered iID as index of the point,  jID the index of the bond of iID and nID the neighborID.
# Returns
- `stress_NP1::SubArray`: updated stresses
Example:
```julia
```
"""

function _elastic_stresses!(nodes::AbstractVector{Int64},
                            dof::Int64,
                            strain_increment::NodeTensorField{Float64},
                            stress_N::NodeTensorField{Float64},
                            stress_NP1::NodeTensorField{Float64})

    mapping = if dof == 2
        get_mapping(2)::SMatrix{3,2,Int64,6}
    elseif dof == 3
        get_mapping(3)::SMatrix{6,2,Int64,12}
    else
        get_mapping(dof)
    end

    hooke_matrix::NodeTensorField{Float64} = Data_Manager.get_field("Material Gradient")

    @inbounds for iID in nodes
        # Voigt-Notation
        for m in axes(mapping, 1)
            i = mapping[m, 1]
            j = mapping[m, 2]

            sNP1 = stress_N[iID, i, j]

            # Constitutive Relation: σ_ij^(n+1) = σ_ij^n + C_ijkl * Δε_kl
            for k in axes(mapping, 1)
                i_k = mapping[k, 1]
                j_k = mapping[k, 2]
                factor = (i_k == j_k) ? 1.0 : 2.0   # Voigt: γ = 2ε für Schub
                sNP1 += hooke_matrix[iID, m, k] *
                        factor * strain_increment[iID, i_k, j_k]
            end

            stress_NP1[iID, i, j] = sNP1
            if i != j
                stress_NP1[iID, j, i] = sNP1
            end
        end
    end
end

"""
    compute_stresses(nodes, dof, p, material, time, dt, strain_increment, stress_N, stress_NP1)

Computes the stresses of `nodes`.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `dof::Int64`: Degrees of freedom.
- `p`: The model parameters (the model's `@params` struct).
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strain_increment`: Strain increment.
- `stress_N`: Stress of step N.
- `stress_NP1`: Stress of step N+1, updated in place.
"""
function compute_stresses(nodes::AbstractVector{Int64}, dof::Int64,
                          p::CorrespondenceElasticParams, material,
                          time::Float64, dt::Float64,
                          strain_increment::NodeTensorField{Float64},
                          stress_N::NodeTensorField{Float64},
                          stress_NP1::NodeTensorField{Float64})
    return _elastic_stresses!(nodes, dof, strain_increment, stress_N, stress_NP1)
end

"""
    init_model(nodes, p, material)

Initializes the model fields of `nodes`.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `p`: The model parameters (the model's `@params` struct).
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
"""
function init_model(nodes::AbstractVector{Int64}, p::CorrespondenceElasticParams, material)
    dof::Int64 = Data_Manager.get_dof()
    hooke::NodeTensorField{Float64} = Data_Manager.create_constant_node_tensor_field("Material Gradient",
                                                                                     Float64,
                                                                                     Int64((dof *
                                                                                            (dof +
                                                                                             1)) /
                                                                                           2))
    for iID in nodes
        @views hooke[iID, :, :] = hooke_matrix(material, dof, iID)
    end
end

"""
    compute_stresses_ba(nodes, nlist, dof, p, material, time, dt, strain_increment, stress_N, stress_NP1)

Computes the bond-associated stresses of `nodes`.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `nlist`: The neighborhood list.
- `dof::Int64`: Degrees of freedom.
- `p`: The model parameters (the model's `@params` struct).
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strain_increment`: Strain increment.
- `stress_N`: Stress of step N.
- `stress_NP1`: Stress of step N+1, updated in place.
"""
function compute_stresses_ba(nodes, nlist, dof::Int64, p::CorrespondenceElasticParams,
                             material, time::Float64, dt::Float64, strain_increment,
                             stress_N, stress_NP1)
    @views mapping = get_mapping(dof)
    for iID in nodes
        @views hookeMatrix = hooke_matrix(material, dof, iID)
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            @views sNP1 = stress_NP1[iID][jID, :, :]
            @views sInc = strain_increment[iID][jID, :, :]
            @views sN = stress_N[iID][jID, :, :]
            fast_mul!(sNP1, hookeMatrix, sInc, sN, mapping)
        end
    end
end


"""
	fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
