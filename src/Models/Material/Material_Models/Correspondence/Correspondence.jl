# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence
using LinearAlgebra
using StaticArrays
using LoopVectorization
using TimerOutputs
using Rotations
include("../Zero_Energy_Control/global_control.jl")
using .Global_Zero_Energy_Control: compute_control
include("./Bond_Associated_Correspondence.jl")
using .Bond_Associated_Correspondence
using ....Material_Basis: compute_Piola_Kirchhoff_stress!
using .......Helpers: invert, rotate, determinant, smat, matrix_diff!, fast_mul!, mat_mul!
using .......Geometry: compute_strain
using .Global_Zero_Energy_Control
include("../../../../Core/Module_inclusion/set_Modules.jl")
global module_list = find_module_files(@__DIR__, "correspondence_name")
include_files(module_list)

export init_model
export material_name
export compute_model
export fields_for_local_synchronization
"""
  init_model(datamanager::Module, nodes::AbstractVector{Int64}, block::Int64, material_parameter::Dict{String,Any})

Initializes the correspondence material model.

# Arguments
  - `datamanager::Data_Manager`: Datamanager.
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `block::Int64`: Block id of the current block.
  - `material_parameter::Dict{String,Any}`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_Manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    block::Int64,
                    material_parameter::Dict{String,Any})
    # global dof
    # global rotation
    # global angles
    if !haskey(material_parameter, "Symmetry")
        @error "Symmetry for correspondence material is missing; options are 'isotropic plane strain', 'isotropic plane stress', 'anisotropic plane stress', 'anisotropic plane stress','isotropic' and 'anisotropic'. For 3D the plane stress or plane strain option is ignored."
        return nothing
    end
    dof = datamanager.get_dof()
    datamanager.create_node_field("Strain", Float64, dof, VectorOrMatrix = "Matrix")
    datamanager.create_constant_node_field("Strain Increment", Float64, dof,
                                           VectorOrMatrix = "Matrix")
    datamanager.create_node_field("Cauchy Stress", Float64, dof, VectorOrMatrix = "Matrix")
    datamanager.create_node_field("von Mises Stress", Float64, 1)
    rotation::Bool = datamanager.get_rotation()
    material_models = split(material_parameter["Material Model"], "+")
    material_models = map(r -> strip(r), material_models)
    #occursin("Correspondence", material_name)
    for material_model in material_models
        datamanager.set_analysis_model("Correspondence Model", block, material_model)
        mod = create_module_specifics(material_model,
                                      module_list,
                                      "correspondence_name")
        if isnothing(mod)
            @error "No correspondence material of name " * material_model * " exists."
            return nothing
        end
        datamanager.set_model_module(material_model, mod)
        datamanager = mod.init_model(datamanager, nodes,
                                     material_parameter)
    end
    if haskey(material_parameter, "Bond Associated") &&
       material_parameter["Bond Associated"]
        return Bond_Associated_Correspondence.init_model(datamanager,
                                                         nodes,
                                                         material_parameter)
    end
    material_parameter["Bond Associated"] = false
    return datamanager
end

"""
    material_name()

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
function material_name()
    return "Correspondence"
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module,
                                          model::String,
                                          block::Int64,
                                          model_param::Dict)
    for material_model in datamanager.get_analysis_model("Correspondence Model", block)
        mod = datamanager.get_model_module(material_model)

        datamanager = mod.fields_for_local_synchronization(datamanager,
                                                           model)
        if model_param["Bond Associated"]
            Bond_Associated_Correspondence.fields_for_local_synchronization(datamanager,
                                                                            model)
        end
    end
    return datamanager
end

"""
    compute_model(datamanager, nodes, material_parameter, time, dt, to::TimerOutput)

Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `material_parameter::Dict{String,Any}`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_Manager`: Datamanager.
Example:
```julia
```
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       material_parameter::Dict{String,Any},
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    if material_parameter["Bond Associated"]
        return Bond_Associated_Correspondence.compute_model(datamanager,
                                                            nodes,
                                                            material_parameter,
                                                            block,
                                                            time,
                                                            dt,
                                                            to)
    end
    return compute_correspondence_model(datamanager,
                                        nodes,
                                        material_parameter,
                                        block,
                                        time,
                                        dt,
                                        to)
end

function compute_correspondence_model(datamanager::Module,
                                      nodes::AbstractVector{Int64},
                                      material_parameter::Dict{String,Any},
                                      block::Int64,
                                      time::Float64,
                                      dt::Float64,
                                      to::TimerOutput)
    rotation::Bool = datamanager.get_rotation()
    dof = datamanager.get_dof()
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    bond_force = datamanager.get_field("Bond Forces")
    bond_damage = datamanager.get_bond_damage("NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")

    strain_N = datamanager.get_field("Strain", "N")
    strain_NP1 = datamanager.get_field("Strain", "NP1")
    stress_N = datamanager.get_field("Cauchy Stress", "N")
    stress_NP1 = datamanager.get_field("Cauchy Stress", "NP1")
    strain_increment = datamanager.get_field("Strain Increment")
    compute_strain(nodes, deformation_gradient, strain_NP1)
    matrix_diff!(strain_increment, nodes, strain_NP1, strain_N)

    if rotation
        rotation_tensor = datamanager.get_field("Rotation Tensor")
        stress_N = rotate(nodes, stress_N, rotation_tensor, false)
        strain_increment = rotate(nodes, strain_increment, rotation_tensor, false)
    end

    # material_models = split(material_parameter["Material Model"], "+")
    # material_models = map(r -> strip(r), material_models)
    @timeit to "Calculated material" begin
        for material_model in datamanager.get_analysis_model("Correspondence Model", block)
            mod = datamanager.get_model_module(material_model)

            mod.compute_stresses(datamanager,
                                 nodes,
                                 dof,
                                 material_parameter,
                                 time,
                                 dt,
                                 strain_increment,
                                 stress_N,
                                 stress_NP1)
        end
    end

    if rotation
        stress_NP1 = rotate(nodes, stress_NP1, rotation_tensor, true)
    end
    @timeit to "Calculated bond force" calculate_bond_force!(nodes,
                                                             dof,
                                                             deformation_gradient,
                                                             undeformed_bond,
                                                             bond_damage,
                                                             inverse_shape_tensor,
                                                             stress_NP1,
                                                             bond_force)
    # TODO general interface, because it might be a flexbile Set_modules interface in future
    @timeit to "zero energy" datamanager=zero_energy_mode_compensation(datamanager, nodes,
                                                                       material_parameter,
                                                                       time, dt)
    return datamanager
end

"""
    zero_energy_mode_compensation(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict{String,Any}, time::Float64, dt::Float64)

Global - J. Wan et al., "Improved method for zero-energy mode suppression in peridynamic correspondence model in Acta Mechanica Sinica https://doi.org/10.1007/s10409-019-00873-y

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `material_parameter::Dict{String, Any}`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function zero_energy_mode_compensation(datamanager::Module,
                                       nodes::AbstractVector{Int64},
                                       material_parameter::Dict{String,Any},
                                       time::Float64,
                                       dt::Float64)
    if !haskey(material_parameter, "Zero Energy Control")
        if time == 0
            @warn "No zero energy control activated for corresponcence."
        end
        return datamanager::Module
    end
    if material_parameter["Zero Energy Control"]::String ==
       Global_Zero_Energy_Control.control_name()::String
        datamanager = Global_Zero_Energy_Control.compute_control(datamanager,
                                                                 nodes,
                                                                 material_parameter,
                                                                 time,
                                                                 dt)::Module
    end
    return datamanager
end

"""
    calculate_bond_force(nodes::AbstractVector{Int64}, deformation_gradient::Array{Float64, 3}, undeformed_bond::Vector{Matrix{Float64}}, bond_damage::Vector{Vector{Float64}}, inverse_shape_tensor::Array{Float64, 3}, stress_NP1::Array{Float64, 3}, bond_force::Vector{Matrix{Float64}})

Calculate bond forces for specified nodes based on deformation gradients.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `deformation_gradient::SubArray`: Deformation gradient.
- `undeformed_bond::Vector{Vector{Vector{Float64}}}`: Undeformed bond geometry.
- `bond_damage::Vector{Vector{Float64}}`: Bond damage.
- `inverse_shape_tensor::Array{Float64, 3}`: Inverse shape tensor.
- `stress_NP1::Array{Float64, 3}`: Stress at time step n+1.
- `bond_force::Vector{Vector{Vector{Float64}}}`: Bond force.
# Returns
- `bond_force::Vector{Vector{Vector{Float64}}}`: Bond force.
"""
function calculate_bond_force!(nodes::AbstractVector{Int64},
                               dof::Int64,
                               deformation_gradient::Array{Float64,3},
                               undeformed_bond::Vector{Vector{Vector{Float64}}},
                               bond_damage::Vector{Vector{Float64}},
                               inverse_shape_tensor::Array{Float64,3},
                               stress_NP1::Array{Float64,3},
                               bond_force::Vector{Vector{Vector{Float64}}})
    if dof == 2
        return calculate_bond_force_2d!(bond_force,
                                        nodes,
                                        deformation_gradient,
                                        undeformed_bond,
                                        bond_damage,
                                        inverse_shape_tensor,
                                        stress_NP1)
    elseif dof == 3
        return calculate_bond_force_3d!(bond_force,
                                        nodes,
                                        deformation_gradient,
                                        undeformed_bond,
                                        bond_damage,
                                        inverse_shape_tensor,
                                        stress_NP1)
    end
    #return bond_force
end

function calculate_bond_force_2d!(bond_force::Vector{Vector{Vector{Float64}}},
                                  nodes::AbstractVector{Int64},
                                  deformation_gradient::Array{Float64,3},
                                  undeformed_bond::Vector{Vector{Vector{Float64}}},
                                  bond_damage::Vector{Vector{Float64}},
                                  inverse_shape_tensor::Array{Float64,3},
                                  stress_NP1::Array{Float64,3})
    pk_stress = MMatrix{2,2}(zeros(2, 2))
    temp = MMatrix{2,2}(zeros(Float64, 2, 2))
    @inbounds @fastmath for iID in nodes
        @views compute_Piola_Kirchhoff_stress!(SMatrix{2,2}(@views stress_NP1[iID, :,
                                                                              :]),
                                               SMatrix{2,2}(@views deformation_gradient[iID,
                                                                                        :,
                                                                                        :]),
                                               pk_stress)

        mat_mul!(temp, SMatrix{2,2}(pk_stress),
                 SMatrix{2,2}(@views inverse_shape_tensor[iID, :, :]))

        @views @inbounds @fastmath for jID in axes(bond_damage[iID], 1)
            @inbounds @fastmath for idof in 1:2
                b_fi = Float64(0.0)
                @inbounds @fastmath for jdof in 1:2
                    b_fi += temp[idof, jdof] *
                            bond_damage[iID][jID] *
                            undeformed_bond[iID][jID][jdof]
                end
                bond_force[iID][jID][idof] = b_fi
            end
        end
    end
end

function calculate_bond_force_3d!(bond_force::Vector{Vector{Vector{Float64}}},
                                  nodes::AbstractVector{Int64},
                                  deformation_gradient::Array{Float64,3},
                                  undeformed_bond::Vector{Vector{Vector{Float64}}},
                                  bond_damage::Vector{Vector{Float64}},
                                  inverse_shape_tensor::Array{Float64,3},
                                  stress_NP1::Array{Float64,3})
    pk_stress = MMatrix{3,3}(zeros(3, 3))
    temp = MMatrix{3,3}(zeros(Float64, 3, 3))
    @inbounds @fastmath for iID in nodes
        @views compute_Piola_Kirchhoff_stress!(SMatrix{3,3}(stress_NP1[iID, :,
                                                                       :]),
                                               SMatrix{3,3}(deformation_gradient[iID,
                                                                                 :,
                                                                                 :]),
                                               pk_stress)

        mat_mul!(temp, SMatrix{3,3}(pk_stress),
                 SMatrix{3,3}(inverse_shape_tensor[iID, :, :]))
        @views @inbounds @fastmath for jID in axes(bond_damage[iID], 1)
            @inbounds @fastmath for idof in 1:3
                b_fi = Float64(0.0)
                @inbounds @fastmath for jdof in 1:3
                    @views b_fi += temp[idof, jdof] *
                                   bond_damage[iID][jID] *
                                   undeformed_bond[iID][jID][jdof]
                end
                bond_force[iID][jID][idof] = b_fi
            end
        end
    end
end

end
