# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Plastic
using ......Data_Manager
using .....Material_Basis:
                           flaw_function, get_von_mises_yield_stress,
                           compute_deviatoric_and_spherical_stresses
using LinearAlgebra
using StaticArrays
export fe_support
export init_model
export fields_for_local_synchronization

yield_stress::Float64 = 1.0
reduced_yield_stress::Float64 = 0.0
spherical_stress_N::Float64 = 0.0
spherical_stress_NP1::Float64 = 0.0

deviatoric_stress_N::Matrix{Float64} = zeros(2, 2)
deviatoric_stress_NP1::Matrix{Float64} = zeros(2, 2)
temp_A::Matrix{Float64} = zeros(2, 2)
temp_B::Matrix{Float64} = zeros(2, 2)

temp_scalar::Float64 = 0.0

sqrt23::Float64 = sqrt(2 / 3)
deviatoric_stress_magnitude_N::Float64 = 0.0
deviatoric_stress_magnitude_NP1::Float64 = 0.0

dev_strain_inc::Matrix{Float64} = zeros(2, 2)
spherical_strain::Float64 = 0.0

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
    return false
end

"""
  init_model(nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
"""
function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    if !haskey(material_parameter, "Shear Modulus")
        @error "Shear Modulus must be defined to be able to run this plastic material"
        return nothing
    end
    if !haskey(material_parameter, "Yield Stress")
        @error "No ''Yield Stress'' is defined."
        return nothing
    end

    Data_Manager.create_node_scalar_field("von Mises Yield Stress", Float64)
    Data_Manager.create_node_scalar_field("Plastic Strain", Float64)

    if haskey(material_parameter, "Bond Associated") &&
       material_parameter["Bond Associated"]
        Data_Manager.create_bond_scalar_state("von Mises Bond Yield Stress", Float64)
        Data_Manager.create_bond_scalar_state("Plastic Bond Strain", Float64)
    end
    dof = Data_Manager.get_dof()
    deviatoric_stress_N = zeros(dof, dof)
    deviatoric_stress_NP1 = zeros(dof, dof)
    temp_A = zeros(dof, dof)
    temp_B = zeros(dof, dof)
    dev_strain_inc = zeros(dof, dof)
end

"""
   correspondence_name()

   Gives the material name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: The name of the material.

   Example:
   ```julia
   println(material_name())
   "Material Template"
   ```
   """
function correspondence_name()
    return "Correspondence Plastic"
end

"""
	compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray, iID_jID_nID::Tuple=())

Calculates the stresses of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `iID::Int64`: Node ID.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{NodeTensorField{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
- `iID_jID_nID::Tuple=(): (optional) are the index and node id information. The tuple is ordered iID as index of the point,  jID the index of the bond of iID and nID the neighborID.
# Returns
- `stress_NP1::SubArray`: updated stresses

Example:
```julia
```
"""
function compute_stresses(nodes,
                          dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::Union{SubArray,NodeTensorField{Float64,3}},
                          stress_N::Union{SubArray,NodeTensorField{Float64,3}},
                          stress_NP1::Union{SubArray,NodeTensorField{Float64,3}})
    global yield_stress
    global reduced_yield_stress
    global spherical_stress_N
    global spherical_stress_NP1
    global deviatoric_stress_N
    global deviatoric_stress_NP1
    global temp_A
    global temp_B
    global sqrt23

    von_Mises_stress_yield::NodeScalarField{Float64} = Data_Manager.get_field("von Mises Yield Stress",
                                                                              "NP1")
    plastic_strain_N::NodeScalarField{Float64} = Data_Manager.get_field("Plastic Strain",
                                                                        "N")
    plastic_strain_NP1::NodeScalarField{Float64} = Data_Manager.get_field("Plastic Strain",
                                                                          "NP1")
    coordinates::NodeVectorField{Float64} = Data_Manager.get_field("Coordinates")
    yield_stress = material_parameter["Yield Stress"]
    spherical_stress_N = 0.0
    spherical_stress_NP1 = 0.0
    deviatoric_stress_N .= 0.0
    deviatoric_stress_NP1 .= 0.0
    temp_A .= 0.0
    temp_B .= 0.0

    # sqrt23::Float64 = sqrt(2 / 3)
    for iID in nodes
        # @views reduced_yield_stress = yield_stress
        reduced_yield_stress = flaw_function(material_parameter, coordinates[iID, :],
                                             yield_stress)

        stress_NP1[iID, :, :], plastic_strain_NP1[iID],
        von_Mises_stress_yield[iID] = compute_plastic_model(stress_NP1[iID,
                                                            :,
                                                            :],
                                                            stress_N[iID,
                                                            :,
                                                            :],
                                                            spherical_stress_NP1,
                                                            spherical_stress_N,
                                                            deviatoric_stress_NP1,
                                                            deviatoric_stress_N,
                                                            strain_increment[iID,
                                                            :,
                                                            :],
                                                            von_Mises_stress_yield[iID],
                                                            plastic_strain_NP1[iID],
                                                            plastic_strain_N[iID],
                                                            reduced_yield_stress,
                                                            material_parameter["Shear Modulus"],
                                                            dof,
                                                            temp_A,
                                                            temp_B,
                                                            sqrt23)
    end
end

function compute_stresses_ba(nodes,
                             nlist,
                             dof::Int64,
                             material_parameter::Dict,
                             time::Float64,
                             dt::Float64,
                             strain_increment::Vector{AbstractArray{Float64,3}},
                             stress_N::Vector{AbstractArray{Float64,3}},
                             stress_NP1::Vector{AbstractArray{Float64,3}})
    temp_A = @MMatrix zeros(dof, dof)
    temp_B = @MMatrix zeros(dof, dof)

    sqrt23::Float64 = sqrt(2 / 3)
    von_Mises_stress_yield = Data_Manager.get_field("von Mises Bond Yield Stress", "NP1")
    plastic_strain_N = Data_Manager.get_field("Plastic Bond Strain", "N")
    plastic_strain_NP1 = Data_Manager.get_field("Plastic Bond Strain", "NP1")
    coordinates = Data_Manager.get_field("Coordinates")
    yield_stress::Float64 = material_parameter["Yield Stress"]
    spherical_stress_N::Float64 = 0
    deviatoric_stress_N = @MMatrix zeros(dof, dof)

    spherical_stress_NP1::Float64 = 0
    deviatoric_stress_NP1 = @MMatrix zeros(dof, dof)

    for iID in nodes
        @views reduced_yield_stress = yield_stress
        @views reduced_yield_stress = flaw_function(material_parameter, coordinates[iID, :],
                                                    yield_stress)
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            stress_NP1[iID][jID, :, :],
            plastic_strain_NP1[iID][jID],
            von_Mises_stress_yield[iID][jID] = compute_plastic_model(stress_NP1[iID][jID, :,
                                                                     :],
                                                                     stress_N[iID][jID, :,
                                                                     :],
                                                                     spherical_stress_NP1,
                                                                     spherical_stress_N,
                                                                     deviatoric_stress_NP1,
                                                                     deviatoric_stress_N,
                                                                     strain_increment[iID][jID,
                                                                     :,
                                                                     :],
                                                                     von_Mises_stress_yield[iID][jID],
                                                                     plastic_strain_NP1[iID][jID],
                                                                     plastic_strain_N[iID][jID],
                                                                     reduced_yield_stress,
                                                                     material_parameter["Shear Modulus"],
                                                                     dof,
                                                                     temp_A,
                                                                     temp_B,
                                                                     sqrt23)
        end
    end
end

function compute_plastic_model(stress_NP1,
                               stress_N,
                               spherical_stress_NP1,
                               spherical_stress_N,
                               deviatoric_stress_NP1,
                               deviatoric_stress_N,
                               strain_increment,
                               von_Mises_stress_yield,
                               plastic_strain_NP1,
                               plastic_strain_N,
                               reduced_yield_stress,
                               shear_modulus,
                               dof,
                               temp_A,
                               temp_B,
                               sqrt23)
    global temp_scalar
    global sqrt23
    global deviatoric_stress_magnitude_N
    global deviatoric_stress_magnitude_NP1
    global dev_strain_inc
    global spherical_strain

    compute_deviatoric_and_spherical_stresses(stress_NP1,
                                              spherical_stress_NP1,
                                              deviatoric_stress_NP1,
                                              dof)

    von_Mises_stress_yield = get_von_mises_yield_stress(deviatoric_stress_NP1)
    if von_Mises_stress_yield < reduced_yield_stress
        # material is elastic and nothing happens
        plastic_strain_NP1 = plastic_strain_N
        return stress_NP1, plastic_strain_NP1, von_Mises_stress_yield
    end
    deviatoric_stress_magnitude_NP1 = maximum([1.0e-20, von_Mises_stress_yield / sqrt23])
    deviatoric_stress_NP1 .*= sqrt23 * reduced_yield_stress /
                              deviatoric_stress_magnitude_NP1
    stress_NP1 = deviatoric_stress_NP1 + spherical_stress_NP1 .* I(dof)

    von_Mises_stress_yield = get_von_mises_yield_stress(deviatoric_stress_NP1)
    #https://de.wikipedia.org/wiki/Plastizit%C3%A4tstheorie
    #############################
    # comment taken from Peridigm elastic_plastic_correspondence.cxx
    #############################
    # Update the equivalent plastic strain
    #
    # The algorithm below is generic and should not need to be modified for
    # any J2 plasticity yield surface.  It uses the difference in the yield
    # surface location at the NP1 and N steps to increment eqps regardless
    # of how the plastic multiplier was found in the yield surface
    # evaluation.
    #
    # First go back to step N and compute deviatoric stress and its
    # magnitude.  We didn't do this earlier because it wouldn't be necassary
    # if the step is elastic.
    compute_deviatoric_and_spherical_stresses(stress_N,
                                              spherical_stress_N,
                                              deviatoric_stress_N,
                                              dof)
    deviatoric_stress_magnitude_N = maximum([1.0e-20, norm(deviatoric_stress_N)])
    # Contract the two tensors. This represents a projection of the plastic
    # strain increment tensor onto the "direction" of deviatoric stress
    # increment
    # deviatoricStrainInc[i]
    dev_strain_inc .= 0
    spherical_strain = 0
    compute_deviatoric_and_spherical_stresses(strain_increment,
                                              spherical_strain,
                                              dev_strain_inc,
                                              dof)

    temp_A = dev_strain_inc -
             (deviatoric_stress_NP1 - deviatoric_stress_N) ./ 2 / shear_modulus
    temp_B = (deviatoric_stress_NP1 ./ deviatoric_stress_magnitude_NP1 +
              deviatoric_stress_N ./ deviatoric_stress_magnitude_N) ./ 2
    temp_scalar = sum(temp_A .* temp_B) # checken

    plastic_strain_NP1 = plastic_strain_N + maximum([0, sqrt23 * temp_scalar])
    return stress_NP1, plastic_strain_NP1, von_Mises_stress_yield
end

"""
	fields_for_local_synchronization( model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
