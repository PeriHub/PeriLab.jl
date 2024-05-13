# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using Statistics
include("../Physics/Material/material_basis.jl")

function get_FE_material_model(params::Dict, name::String)
    if !haskey(params["Material Models"], params["FEM"][name]["Material Model"])
        @error "Material model " * params["FEM"][name]["Material Model"] * " defined in FEM are not defined as material"
        return nothing
    end
    return params["Material Models"][params["FEM"][name]["Material Model"]]
end

function calculate_FEM(datamanager::Module, elements::Union{SubArray,Vector{Int64}}, params::Dict, compute_stresses, time::Float64, dt::Float64)

    rotation::Bool, angles = datamanager.rotation_data("Element")
    dof = datamanager.get_dof()

    force_densities = datamanager.get_field("Force Densities", "NP1")
    displacement = datamanager.get_field("Displacements", "NP1")
    strain_N = datamanager.get_field("Element Strain", "N")
    strain_NP1 = datamanager.get_field("Element Strain", "NP1")
    stress_N = datamanager.get_field("Element Stress", "N")
    stress_NP1 = datamanager.get_field("Element Stress", "NP1")
    strain_increment = datamanager.get_field("Element Strain Increment")
    topology = datamanager.get_field("FE Topology")
    jacobian = datamanager.get_field("Element Jacobi Matrix")
    det_jacobian = datamanager.get_field("Element Jacobi Determinant")
    # not avveraged Cauchy stresses -> element stresses will be later available in Exodus
    cauchy_stress = datamanager.get_field("Cauchy Stress", "NP1")

    B_matrix = datamanager.get_field("B Matrix")

    stress_temp = @MVector zeros(3 * dof - 3)
    le::Int64 = 0
    for id_el in elements
        topo = view(topology, id_el, :)
        le = dof * length(topo)
        nnodes::Int64 = length(topo)

        for id_int in eachindex(B_matrix[:, 1, 1])

            strain_NP1[id_el, id_int, :] = B_matrix[id_el, id_int, :, :]' * reshape((displacement[topo, :])', le)
            strain_increment[id_el, id_int, :] = strain_NP1[id_el, id_int, :] - strain_N[id_el, id_int, :]

            if rotation
                #tbd
                stress_N = rotate(nodes, dof, stress_N, angles, false)
                strain_increment = rotate(nodes, dof, strain_increment, angles, false)
            end

            # in future this part must be changed -> using set Modules

            stress_NP1[id_el, id_int, :], datamanager = compute_stresses(datamanager, dof, params["Material Model"], time, dt, strain_increment[id_el, id_int, :], stress_N[id_el, id_int, :], stress_NP1[id_el, id_int, :])

            #specifics = Dict{String,String}("Call Function" => "compute_stresses", "Name" => "material_name") -> tbd
            # material_model is missing
            #stress_NP1, datamanager = Set_modules.create_module_specifics(material_model, module_list, specifics, (datamanager, nodes, dof, material_parameter, time, dt, strain_increment, stress_N, stress_NP1))

            if rotation
                #tbd
                stress_NP1 = rotate(nodes, dof, stress_NP1, angles, true)
            end

            force_densities[topo, :] += reshape(B_matrix[id_el, id_int, :, :] * stress_NP1[id_el, id_int, :] .* det_jacobian[id_el, id_int], (dof, nnodes))'
            # if you do not use permutedims you will get some index errors
            stress_temp .+= stress_NP1[id_el, id_int, :] .* det_jacobian[id_el, id_int]
        end
        # as long as no elements stresses are written
        temp = permutedims(cauchy_stress[topo, :, :], (2, 3, 1))
        temp[:, :, 1:nnodes] .= voigt_to_matrix(stress_temp)
        cauchy_stress[topo, :, :] = permutedims(temp[:, :, 1:nnodes], (3, 1, 2)) ./ sum(det_jacobian[id_el, :])
        stress_temp .= 0
    end
    return datamanager
end

function get_lumped_mass(elements::Vector{Int64}, dof::Int64, topology::SubArray{Int64}, N::SubArray{Float64}, determinant_jacobian::SubArray{Float64}, rho::SubArray{Float64}, lumped_mass::SubArray{Float64})
    for id_int in eachindex(N[:, 1, 1])
        temp = N[id_int, :, :] * N[id_int, :, :]'
        for id_el in elements
            nnodes = length(topology[id_el, :])
            mean_rho = mean(rho[topology[id_el, :]])
            for i_node in 1:nnodes
                lumped_mass[topology[id_el, i_node]] += sum(temp[(i_node-1)*dof+1, :]) .* mean_rho * determinant_jacobian[id_el, id_int]
            end
        end
    end
    return lumped_mass
end



"""
    get_weights_and_integration_points(dof::Int64, p::Vector{Int64})

Compute integration points and weights using Gauss-Legendre quadrature.

# Arguments
- `dof::Int64`: The number of degrees of freedom.
- `p::Vector{Int64}`: A vector containing the polynomial degrees for each degree of freedom.

# Returns
- `w::Matrix{Float64}`: Matrix of integration weights. Each row corresponds to a degree of freedom, and columns contain the integration weights.
- `x::Matrix{Float64}`: Matrix of integration points. Each row corresponds to a degree of freedom, and columns contain the integration points.

# Example
```julia
dof = 3
p = [2, 3, 1]
x, w = get_weights_and_integration_points(dof, p)
"""
function get_weights_and_integration_points(dof::Int64, num_int::Vector{Int64})

    x::Matrix{Float64} = zeros(Float64, maximum(num_int), dof)
    w::Matrix{Float64} = zeros(Float64, maximum(num_int), dof)
    for idof in 1:dof
        x[1:num_int[idof], idof], w[1:num_int[idof], idof] = gausslegendre(num_int[idof])
    end
    return w, x
end

"""
    get_multi_dimensional_integration_point_data(dof::Int64, p::Vector{Int64}, value::Matrix{Float64})

Restructure integration point information for multi-dimensional problems.

# Arguments
- `dof::Int64`: Degree of freedom, only 2 and 3 are supported.
- `p::Vector{Int64}`: Vector containing the number of integration points in each dimension.
- `value::Matrix{Float64}`: Matrix containing the coordinates of the reference element.

# Returns
- `integration_point::Matrix{Float64}`: Matrix of integration point data.

# Example
```julia
dof = 2
p = [2, 3]
value = rand(3, dof)
result = get_multi_dimensional_integration_points(dof, p, value)
"""
function get_multi_dimensional_integration_point_data(dof::Int64, num_int::Vector{Int64}, value::Matrix{Float64})
    count::Int64 = 0
    integration_point::Matrix{Float64} = zeros(prod(num_int), dof)

    if dof == 2
        for jID in 1:num_int[2]
            for iID in 1:num_int[1]
                count += 1
                integration_point[count, 1] = value[iID, 1]
                integration_point[count, 2] = value[jID, 2]
            end
        end
    elseif dof == 3
        for kID in 1:num_int[3]
            for jID in 1:num_int[2]
                for iID in 1:num_int[1]
                    count += 1
                    integration_point[count, 1] = value[iID, 1]
                    integration_point[count, 2] = value[jID, 2]
                    integration_point[count, 3] = value[kID, 3]
                end
            end
        end
    else
        @error "degree of freedom = $dof is not supported, only 2 and 3."
        return nothing
    end
    return integration_point
end


function get_Jacobian(elements::Vector{Int64}, dof::Int64, topology::SubArray{Int64}, coordinates::Union{SubArray{Float64},SubArray{Float64}}, B::Union{Array{Float64},SubArray}, jacobian::SubArray{Float64}, determinant_jacobian::SubArray{Float64})

    mapping = Vector{Int64}(1:dof:length(B[1, :, 1]))
    for id_el in elements
        for id_int in eachindex(B[:, 1, 1])
            for idof in 1:dof
                for jdof in 1:dof

                    jacobian[id_el, id_int, idof, jdof] = dot(coordinates[topology[id_el, :], idof], B[id_int, mapping.+(jdof-1), jdof])
                end
            end

            determinant_jacobian[id_el, id_int] = det(jacobian[id_el, id_int, :, :])
            if determinant_jacobian[id_el, id_int] <= 0
                @error "The determinat of the Jacobian is " * string(determinant_jacobian[id_el, id_int]) * " in local element $id_el, and must be greater zero."
                return nothing, nothing
            end
            jacobian[id_el, id_int, :, :] = inv(jacobian[id_el, id_int, :, :])
        end
    end
    return jacobian, determinant_jacobian
end


function get_polynomial_degree(params::Dict, dof::Int64)
    if !haskey(params, "Degree")
        @error "No element degree defined"
        return nothing
    end
    value = params["Degree"]
    if typeof(value) == String
        value = parse.(Float64, split(value))
    end
    if sum(typeof.(value) .!= Int64) != 0
        value = Int64.(round.(value))
    end
    if length(value) == 1
        return_value::Vector{Int64} = zeros(dof)
        return_value[1:dof] .= value[1]
        return return_value
    elseif length(value) == dof
        return value[1:dof]
    else
        @error "Degree must be defined with length one or number of dof."
        return nothing
    end
end


"""
    create_element_matrices()

Compute the matrix of shape functions (N) and its derivative (B) for 2D and 3D.
N^TN*\rho give than the mass matrix and B^TCB the stiffness matrix [WillbergC2013](@cite)

# Arguments
- `dof::Int64`: The number of degrees of freedom.
- `p::Vector{Int64}`: A vector containing the polynomial degrees for each degree of freedom.

# Returns
- N = [N_x 0   0   ...
       0   N_y 0   ...
       0   0   N_z ...]

- B = [B_x 0   0   ...
       0   B_y 0   ...
       0   0   B_z ...]
       B_y B_x 0   ...]
       B_z 0 B_x   ...]
       0   B_z B_y ...]
# Example

"""

function create_element_matrices(dof::Int64, p::Vector{Int64}, create_matrices)
    num_int = get_number_of_integration_points(p, dof)
    weights, integration_points = get_weights_and_integration_points(dof, num_int)
    ip_coordinates = get_multi_dimensional_integration_point_data(dof, num_int, integration_points)

    if isnothing(ip_coordinates)
        return nothing, nothing
    end
    ip_weights = get_multi_dimensional_integration_point_data(dof, num_int, weights)
    return create_matrices(dof, num_int, p, ip_weights, ip_coordinates)
end



function get_number_of_integration_points(p::Vector{Int64}, dof::Int64)
    num_int::Vector{Int64} = zeros(Int64, dof)
    for idof in 1:dof
        if p[idof] == 1
            num_int[idof] = 2
            continue
        end
        num_int[idof] = 2 * p[idof] - 1
    end
    return num_int
end


function create_B_matrix(elements::Vector{Int64}, dof::Int64, B_elem, jacobian, B)
    # B is (num_elem, num_integration_points, dof*nodes, 3*dof-2)
    int_points = length(B[1, :, 1, 1])
    nnodes = length(B[1, 1, :, 1]) / dof
    for id_el in elements
        for id_int in 1:int_points
            for id_nodes in 1:nnodes
                B[id_el, id_int, (id_nodes-1)*dof+1:(id_nodes)*dof, :] = jacobian[id_el, id_int, :, :] * B[point_id, (id_nodes-1)*dof:(id_nodes)*dof, :]
            end
        end
    end
    return B
end