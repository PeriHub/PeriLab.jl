# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module FEM_Basis

using TimerOutputs: @timeit
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using Statistics

using ....Data_Manager
using ....PeriLabExceptions: @abort
using ....Helpers: invert, determinant, voigt_to_matrix!

function get_FE_material_model(params::Dict{String,Any}, name::String)
    if !haskey(params["Material Models"], params["FEM"][name]["Material Model"])
        @abort "Material model " *
               params["FEM"][name]["Material Model"] *
               " defined in FEM are not defined as material"
        return
    end
    parameter = params["Material Models"][params["FEM"][name]["Material Model"]]

    return parameter
end

function compute_FEM(elements::AbstractVector{Int64},
                     params::Dict{String,Any},
                     compute_stresses!::Function,
                     time::Float64,
                     dt::Float64)
    @timeit "load_data" begin
        dof = Data_Manager.get_dof()

        forces::NodeVectorField{Float64} = Data_Manager.get_field("Forces", "NP1")
        uNP1::NodeVectorField{Float64} = Data_Manager.get_field("Displacements", "NP1")
        strain_N::FreeSizeField{Float64} = Data_Manager.get_field("Element Strain", "N")
        strain_NP1::FreeSizeField{Float64} = Data_Manager.get_field("Element Strain", "NP1")
        stress_N::FreeSizeField{Float64} = Data_Manager.get_field("Element Stress", "N")
        stress_NP1::FreeSizeField{Float64} = Data_Manager.get_field("Element Stress", "NP1")
        strain_increment::FreeSizeField{Float64} = Data_Manager.get_field("Element Strain Increment")
        topology::FreeSizeField{Int64} = Data_Manager.get_field("FE Topology")
        jacobian::FreeSizeField{Float64} = Data_Manager.get_field("Element Jacobi Matrix")
        det_jacobian::FreeSizeField{Float64} = Data_Manager.get_field("Element Jacobi Determinant")
        cauchy_stress::NodeTensorField{Float64} = Data_Manager.get_field("Cauchy Stress",
                                                                         "NP1")
        hooke_matrix::NodeTensorField{Float64} = Data_Manager.get_field("Material Gradient")
        B_matrix::FreeSizeField{Float64} = Data_Manager.get_field("B Matrix")
        connected_el::NodeScalarField{Int64} = Data_Manager.get_field("Connected Elements")
        stress_temp = @MVector zeros(3 * dof - 3)
        stress_temp_matrix = @MMatrix zeros(dof, dof)
        le::Int64 = 0
        hm = zeros(3 * dof - 3, 3 * dof - 3)
    end
    @timeit "eval loop" begin
        for id_el in elements
            topo = view(topology, id_el, :)
            le = dof * length(topo)

            @views hm = avg_mat(hooke_matrix[topo, :, :])
            f_workspace .= zeros(le)

            for id_int in eachindex(B_matrix[1, :, 1, 1])
                @timeit "strain" begin
                    sNP1 = view_function(strain_NP1, id_el, id_int)
                    compute_strain!(sNP1, B_matrix, uNP1, topo, id_el, id_int)
                    sInc = view_function(strain_increment, id_el, id_int)
                    sN = view_function(strain_N, id_el, id_int)
                    diff_strain!(sInc, sNP1, sN)
                end
                stressN = view_function(strain_N, id_el, id_int)
                stressNP1 = view_function(strain_NP1, id_el, id_int)
                @timeit "compute_stresses" compute_stresses!(dof,
                                                             hm,
                                                             time,
                                                             dt,
                                                             sInc,
                                                             stressN,
                                                             stressNP1)

                @timeit "accumulate_forces" accumulate_forces!(forces, f_workspace,
                                                               @view(B_matrix[id_el, id_int,
                                                                              :, :]),
                                                               stressNP1, topo, dof,
                                                               det_jacobian[id_el, id_int])

                @timeit "stress_temp" @views stress_temp .+= stress_NP1[id_el, id_int, :] .*
                                                             det_jacobian[id_el, id_int]
            end

            # as long as no elements stresses are written
            @timeit "voigt_to_matrix" begin
                voigt_to_matrix!(stress_temp_matrix, stress_temp)
                add_Cauchy_stress!(cauchy_stress, stress_temp_matrix, dof, topo,
                                   connected_el)
            end

            stress_temp .= 0
        end
    end
end

function compute_strain!(strain::AbstractVector{Float64},
                         B::Array{Float64,4},   # (le, nstrain), le = nnodes*ndof
                         uNP1::Matrix{Float64}, # (:, ndof)
                         topo::AbstractVector{<:Integer}, id_el::Int64, id_int::Int64)
    nstrain = size(B, 2)
    nnodes = length(topo)
    ndof = size(uNP1, 2)
    @inbounds for s in 1:nstrain
        acc = 0.0
        for n in 1:nnodes
            node = topo[n]
            base = (n - 1) * ndof
            for d in 1:ndof
                acc += B[id_el, id_int, base + d, s] * uNP1[node, d]
            end
        end
        strain[s] = acc
    end
end

function add_Cauchy_stress!(cauchy_stress::NodeTensorField{Float64},
                            stress_temp_matrix::AbstractMatrix{Float64},
                            dof::Int64,
                            topo::AbstractVector{Int64},
                            connected_el::NodeScalarField{Int64})
    for node in topo
        for i in 1:dof
            for j in 1:dof
                cauchy_stress[node, i, j] += stress_temp_matrix[i, j] / connected_el[node]
            end
        end
    end
end

function accumulate_forces!(forces::AbstractMatrix{Float64},
                            f_workspace::AbstractVector{Float64},  # länge dof*nnodes, vorallokiert
                            B::AbstractMatrix{Float64},
                            stressNP1::AbstractVector{Float64},
                            topo::AbstractVector{<:Integer},
                            dof::Int,
                            scale::Float64)
    mul!(f_workspace, B, stressNP1)
    nnodes = length(topo)
    @inbounds for n in 1:nnodes
        node = topo[n]
        base = (n - 1) * dof
        for d in 1:dof
            forces[node, d] -= f_workspace[base + d] * scale
        end
    end
    return forces
end

function view_function(A::AbstractArray{Float64,3}, i::Int64, j::Int64)
    view(A, i, j, :)
end

function diff_strain!(strain_increment::AbstractVector{Float64},
                      strain_NP1::AbstractVector{Float64},
                      strain_N::AbstractVector{Float64})
    strain_increment .= strain_NP1 .- strain_N
end

function avg_mat(A::AbstractArray{T,3}) where {T}
    # creates the material average matrix for the element, because the material is defined at the points
    dropdims(mean(A, dims = 1), dims = 1)
end

function get_lumped_mass(elements::Vector{Int64},
                         dof::Int64,
                         topology::Matrix{Int64},
                         N::Array{Float64,3},
                         determinant_jacobian::Matrix{Float64},
                         rho::Vector{Float64},
                         lumped_mass::Vector{Float64})
    for id_int in eachindex(N[:, 1, 1])
        temp = N[id_int, :, :] * N[id_int, :, :]'
        for id_el in elements
            nnodes = length(topology[id_el, :])
            mean_rho = mean(rho[topology[id_el, :]])
            for i_node in 1:nnodes
                lumped_mass[topology[id_el, i_node]] += sum(temp[(i_node - 1) * dof + 1, :]) .*
                                                        mean_rho
                # no volume is needed, because the time integration is done F/V
                #* determinant_jacobian[id_el, id_int]
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
```
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
```
"""
function get_multi_dimensional_integration_point_data(dof::Int64,
                                                      num_int::Vector{Int64},
                                                      value::Matrix{Float64})
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
        @abort "degree of freedom = $dof is not supported, only 2 and 3."
        return
    end
    return integration_point
end

function get_Jacobian(elements::Vector{Int64},
                      dof::Int64,
                      topology::Matrix{Int64},
                      coordinates::Matrix{Float64},
                      B::Array{Float64,3},
                      jacobian::Array{Float64,4},
                      determinant_jacobian::Matrix{Float64})
    mapping = Vector{Int64}(1:dof:length(B[1, :, 1]))
    for id_el in elements
        for id_int in eachindex(B[:, 1, 1])
            for idof in 1:dof
                for jdof in 1:dof
                    jacobian[id_el, id_int, idof,
                    jdof] = dot(coordinates[topology[id_el,
                                                              :],
                                                              idof],
                                                              B[id_int,
                                                              mapping .+ (jdof - 1),
                                                              jdof])
                end
            end

            determinant_jacobian[id_el, id_int] = determinant(jacobian[id_el, id_int, :, :])
            if determinant_jacobian[id_el, id_int] <= 0
                @abort "The determinant of the Jacobian is " *
                       string(determinant_jacobian[id_el, id_int]) *
                       " in local element $id_el, and must be greater zero."
                return
            end
            jacobian[id_el, id_int, :,
            :] = invert(jacobian[id_el, id_int, :, :],
                                                   "Jacobian in FEM Module is singular.")
        end
    end
    return jacobian, determinant_jacobian
end

function get_polynomial_degree(params::Dict{String,Any}, dof::Int64)
    if !haskey(params, "Degree")
        @abort "No element degree defined"
        return
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
        @abort "Degree must be defined with length one or number of dof."
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
    ip_coordinates = get_multi_dimensional_integration_point_data(dof, num_int,
                                                                  integration_points)

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

function create_B_matrix(elements::Vector{Int64},
                         dof::Int64,
                         B_elem::Array,
                         jacobian::Array{Float64,4},
                         B::Array{Float64,4})
    # size of B (num_elem, num_integration_points, dof*nodes, 3*dof-2)
    nnodes::Int64 = length(B[1, 1, :, 1]) / dof
    for id_el in elements
        for id_int in eachindex(B[1, :, 1, 1])
            for id_nodes in 1:nnodes
                B[id_el, id_int, ((id_nodes - 1) * dof + 1):((id_nodes) * dof),
                :] = jacobian[id_el,
                                                                                              id_int,
                                                                                              :,
                                                                                              :] *
                                                                                     B_elem[id_int,
                                                                                            ((id_nodes - 1) * dof + 1):(id_nodes * dof),
                                                                                            :]
            end
        end
    end
    return B
end
end
