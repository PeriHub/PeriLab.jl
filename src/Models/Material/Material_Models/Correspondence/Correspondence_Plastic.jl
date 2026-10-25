# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Plastic
using ......Data_Manager
using .......ParameterSpec: @params, Dependent, register_material, value
using ......PeriLabExceptions: @abort
using .....Material_Basis:
                           flaw_function, get_von_mises_yield_stress,
                           compute_deviatoric_and_spherical_stresses
using LinearAlgebra
using StaticArrays
export fe_support
export init_model
export fields_for_local_synchronization

@params "Elastic-plastic correspondence material with von Mises yield criterion" struct CorrespondencePlasticParams
    "Von Mises stress at which the material yields; a number or a table over a node field"
    yield_stress::Dependent = req("Yield Stress"; min = 0, quantity = :stress)
end
__init__() = register_material("Correspondence Plastic", CorrespondencePlasticParams)

const yield_stress::Float64 = 1.0
const reduced_yield_stress::Float64 = 0.0
const spherical_stress_N::Float64 = 0.0
const spherical_stress_NP1::Float64 = 0.0

const deviatoric_stress_N2D::Matrix{Float64} = zeros(2, 2)
const deviatoric_stress_NP12D::Matrix{Float64} = zeros(2, 2)
const dev_strain_inc2D::Matrix{Float64} = zeros(2, 2)

const deviatoric_stress_N3D::Matrix{Float64} = zeros(3, 3)
const deviatoric_stress_NP13D::Matrix{Float64} = zeros(3, 3)
const dev_strain_inc3D::Matrix{Float64} = zeros(3, 3)

const temp_scalar::Float64 = 0.0

const sqrt23::Float64 = sqrt(2 / 3)

const spherical_strain::Float64 = 0.0

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
    init_model(nodes, p, material, block)

Initializes the model fields of `nodes`.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `p`: The model parameters (the model's `@params` struct).
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
"""
function init_model(nodes::AbstractVector{Int64}, p::CorrespondencePlasticParams, material,
                    block::Int64)
    if material.moduli === nothing
        @abort "Shear Modulus must be defined to be able to run this plastic material"
        return
    end
    Data_Manager.create_node_scalar_field("von Mises Yield Stress", Float64)
    Data_Manager.create_node_scalar_field("Plastic Strain", Float64)
    if material.base.bond_associated
        Data_Manager.create_bond_scalar_state("von Mises Bond Yield Stress", Float64)
        Data_Manager.create_bond_scalar_state("Plastic Bond Strain", Float64)
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
function compute_stresses(nodes,
                          dof::Int64,
                          p::CorrespondencePlasticParams,
                          material,
                          time::Float64,
                          dt::Float64,
                          strain_increment::Union{SubArray,NodeTensorField{Float64}},
                          stress_N::Union{SubArray,NodeTensorField{Float64}},
                          stress_NP1::Union{SubArray,NodeTensorField{Float64}})
    if dof == 2
        deviatoric_stress_N = deviatoric_stress_N2D
        deviatoric_stress_NP1 = deviatoric_stress_NP12D
        dev_strain_inc=dev_strain_inc2D
    else
        deviatoric_stress_N = deviatoric_stress_N3D
        deviatoric_stress_NP1 = deviatoric_stress_NP13D
        dev_strain_inc=dev_strain_inc3D
    end

    von_Mises_stress_yield::NodeScalarField{Float64} = Data_Manager.get_field("von Mises Yield Stress",
                                                                              "NP1")
    plastic_strain_N::NodeScalarField{Float64} = Data_Manager.get_field("Plastic Strain",
                                                                        "N")
    plastic_strain_NP1::NodeScalarField{Float64} = Data_Manager.get_field("Plastic Strain",
                                                                          "NP1")
    coordinates::NodeVectorField{Float64} = Data_Manager.get_field("Coordinates")

    # sqrt23::Float64 = sqrt(2 / 3)
    for iID in nodes
        yield_stress = value(p.yield_stress, iID)
        # @views reduced_yield_stress = yield_stress
        reduced_yield_stress = flaw_function(material.base.flaw_function,
                                             view(coordinates, iID, :), yield_stress)
        # the stress of the node is updated in place
        plastic_strain_NP1[iID],
        von_Mises_stress_yield[iID] = compute_plastic_model(view(stress_NP1, iID, :, :),
                                                            view(stress_N, iID, :, :),
                                                            spherical_stress_NP1,
                                                            spherical_stress_N,
                                                            deviatoric_stress_NP1,
                                                            deviatoric_stress_N,
                                                            view(strain_increment, iID, :,
                                                                 :),
                                                            von_Mises_stress_yield[iID],
                                                            plastic_strain_NP1[iID],
                                                            plastic_strain_N[iID],
                                                            reduced_yield_stress,
                                                            material.moduli.shear_modulus,
                                                            dof,
                                                            sqrt23,
                                                            dev_strain_inc)
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
function compute_stresses_ba(nodes,
                             nlist,
                             dof::Int64,
                             p::CorrespondencePlasticParams,
                             material,
                             time::Float64,
                             dt::Float64,
                             strain_increment::Vector{AbstractArray{Float64,3}},
                             stress_N::Vector{AbstractArray{Float64,3}},
                             stress_NP1::Vector{AbstractArray{Float64,3}})
    if dof == 2
        deviatoric_stress_N = deviatoric_stress_N2D
        deviatoric_stress_NP1 = deviatoric_stress_NP12D
        dev_strain_inc=dev_strain_inc2D
    else
        deviatoric_stress_N = deviatoric_stress_N3D
        deviatoric_stress_NP1 = deviatoric_stress_NP13D
        dev_strain_inc=dev_strain_inc3D
    end

    sqrt23::Float64 = sqrt(2 / 3)
    von_Mises_stress_yield = Data_Manager.get_field("von Mises Bond Yield Stress", "NP1")
    plastic_strain_N = Data_Manager.get_field("Plastic Bond Strain", "N")
    plastic_strain_NP1 = Data_Manager.get_field("Plastic Bond Strain", "NP1")
    coordinates = Data_Manager.get_field("Coordinates")
    spherical_stress_N::Float64 = 0
    deviatoric_stress_N = @MMatrix zeros(dof, dof)

    spherical_stress_NP1::Float64 = 0
    deviatoric_stress_NP1 = @MMatrix zeros(dof, dof)

    for iID in nodes
        yield_stress = value(p.yield_stress, iID)
        @views reduced_yield_stress = yield_stress
        @views reduced_yield_stress = flaw_function(material.base.flaw_function,
                                                    coordinates[iID, :],
                                                    yield_stress)
        @fastmath @inbounds for jID in eachindex(nlist[iID])
            # the bond stress is updated in place
            plastic_strain_NP1[iID][jID],
            von_Mises_stress_yield[iID][jID] = compute_plastic_model(view(stress_NP1[iID],
                                                                          jID, :, :),
                                                                     view(stress_N[iID],
                                                                          jID,
                                                                          :, :),
                                                                     spherical_stress_NP1,
                                                                     spherical_stress_N,
                                                                     deviatoric_stress_NP1,
                                                                     deviatoric_stress_N,
                                                                     view(strain_increment[iID],
                                                                          jID, :, :),
                                                                     von_Mises_stress_yield[iID][jID],
                                                                     plastic_strain_NP1[iID][jID],
                                                                     plastic_strain_N[iID][jID],
                                                                     reduced_yield_stress,
                                                                     material.moduli.shear_modulus,
                                                                     dof,
                                                                     sqrt23,
                                                                     dev_strain_inc)
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
                               sqrt23,
                               dev_strain_inc)
    compute_deviatoric_and_spherical_stresses(stress_NP1, spherical_stress_NP1,
                                              deviatoric_stress_NP1, dof)

    von_Mises_stress_yield = get_von_mises_yield_stress(deviatoric_stress_NP1)
    if von_Mises_stress_yield < reduced_yield_stress
        # material is elastic and nothing happens
        plastic_strain_NP1 = plastic_strain_N
        return plastic_strain_NP1, von_Mises_stress_yield
    end
    deviatoric_stress_magnitude_NP1 = max(1.0e-20, von_Mises_stress_yield / sqrt23)
    deviatoric_stress_NP1 .*= sqrt23 * reduced_yield_stress /
                              deviatoric_stress_magnitude_NP1
    # stress_NP1 = deviatoric_stress_NP1 + spherical_stress_NP1 * I, in place
    # (spherical_stress_NP1 is the value passed in, compute_deviatoric_and_spherical_stresses
    # does not return the spherical part)
    stress_NP1 .= deviatoric_stress_NP1
    stress_NP1[diagind(stress_NP1)] .+= spherical_stress_NP1

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
    deviatoric_stress_magnitude_N = max(1.0e-20, norm(deviatoric_stress_N))
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

    # sum(temp_A .* temp_B) with
    # temp_A = dev_strain_inc - (deviatoric_stress_NP1 - deviatoric_stress_N) ./ 2 / shear_modulus
    # temp_B = (deviatoric_stress_NP1 ./ magnitude_NP1 + deviatoric_stress_N ./ magnitude_N) ./ 2
    temp_scalar = 0.0
    for j in axes(dev_strain_inc, 2), i in axes(dev_strain_inc, 1)
        a = dev_strain_inc[i, j] -
            (deviatoric_stress_NP1[i, j] - deviatoric_stress_N[i, j]) / 2 / shear_modulus
        b = (deviatoric_stress_NP1[i, j] / deviatoric_stress_magnitude_NP1 +
             deviatoric_stress_N[i, j] / deviatoric_stress_magnitude_N) / 2
        temp_scalar += a * b
    end

    plastic_strain_NP1 = plastic_strain_N + max(0, sqrt23 * temp_scalar)
    return plastic_strain_NP1, von_Mises_stress_yield
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
