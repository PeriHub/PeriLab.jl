# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Elastic
include("../../Material_Basis.jl")
include("../../../../Support/Helpers.jl")
using .Material_Basis: voigt_to_matrix, matrix_to_voigt, get_Hooke_matrix
using .Helpers: get_fourth_order, fast_mul!, get_mapping
export compute_stresses
export correspondence_name
export fe_support
export init_model
export fields_for_local_synchronization

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
  init_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    return datamanager
end
"""
    correspondence_name()

Gives the material name. It is needed for comparison with the yaml input deck.

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
    compute_stresses(datamanager::Module, iID:Int64, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

Calculates the stresses of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `iID::Int64`: Node ID.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
- `iID_jID_nID::Tuple=(): (optional) are the index and node id information. The tuple is ordered iID as index of the point,  jID the index of the bond of iID and nID the neighborID.
# Returns
- `datamanager::Data_manager`: Datamanager.
- `stress_NP1::SubArray`: updated stresses
Example:
```julia
```
"""
function compute_stresses(datamanager::Module,
                          nodes,
                          dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::AbstractArray{Float64},
                          stress_N::AbstractArray{Float64},
                          stress_NP1::AbstractArray{Float64})
    mapping = get_mapping(dof)
    symmetry = get(material_parameter, "Symmetry", "default")::String
    for iID in nodes
        hookeMatrix = get_Hooke_matrix(datamanager,
                                       material_parameter,
                                       symmetry,
                                       dof,
                                       iID)
        @views sNP1 = stress_NP1[iID, :, :]
        @views sInc = strain_increment[iID, :, :]
        @views sN = stress_N[iID, :, :]
        fast_mul!(sNP1, hookeMatrix, sInc, sN, mapping)
    end
end

function compute_stresses_ba(datamanager::Module,
                             nodes,
                             nlist,
                             dof::Int64,
                             material_parameter::Dict,
                             time::Float64,
                             dt::Float64,
                             strain_increment,
                             stress_N,
                             stress_NP1)
    @views mapping = get_mapping(dof)
    for iID in nodes
        @views hookeMatrix = get_Hooke_matrix(datamanager,
                                              material_parameter,
                                              material_parameter["Symmetry"],
                                              dof,
                                              iID)
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            @views sNP1 = stress_NP1[iID][jID, :, :]
            @views sInc = strain_increment[iID][jID, :, :]
            @views sN = stress_N[iID][jID, :, :]
            fast_mul!(sNP1, hookeMatrix, sInc, sN, mapping)
        end
    end
end

"""
    compute_stresses(datamanager::Module, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

Calculates the stresses of a single node. Needed for FEM. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
# Returns
- `datamanager::Data_manager`: Datamanager.
- `stress_NP1::SubArray`: updated stresses
Example:
```julia
```
"""
function compute_stresses(datamanager::Module,
                          dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::Vector{Float64},
                          stress_N::Vector{Float64},
                          stress_NP1::Vector{Float64})
    hookeMatrix = get_Hooke_matrix(datamanager,
                                   material_parameter,
                                   material_parameter["Symmetry"],
                                   dof)

    return hookeMatrix * strain_increment + stress_N, datamanager
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

end
