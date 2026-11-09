# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

# Test data for PD_Matrix without PeriLab and without the licensed elements: a linear
# spring element (bond-based, own rows only) defined here.

using ArraysOfArrays: VectorOfVectors, flatview, internal_element_ptr
using LinearAlgebra: norm, mul!, Diagonal
using Random
using PeriLab.Solver_Manager.PD_Matrix
import PeriLab.Solver_Manager.PD_Matrix: element_dofs, compute_element!, n_elements

"Spring element: `k_ij = ω_ij V_j / |ξ|³ ξ⊗ξ`, `K[i,i] += k_ij`, `K[i,j] -= k_ij`."
struct TestSpring{D} <: AbstractElement end

n_elements(::TestSpring, data) = length(data.nlist_ptr) - 1

function element_dofs(::TestSpring{D}, data) where {D}
    dofs(i) = collect(((i - 1) * D + 1):(i * D))
    rows = VectorOfVectors{Int64}()
    cols = VectorOfVectors{Int64}()
    for i in 1:(length(data.nlist_ptr) - 1)
        push!(rows, dofs(i))
        push!(cols,
              vcat(dofs(i),
                   [dofs(data.nlist[b])
                    for b in data.nlist_ptr[i]:(data.nlist_ptr[i + 1] - 1)]...))
    end
    return rows, cols
end

@inline function compute_element!(::TestSpring{D}, ke, start, data, i) where {D}
    x = data.coordinates
    first_bond = data.nlist_ptr[i]
    n_bonds = data.nlist_ptr[i + 1] - first_bond
    for q in 0:(D * D * (n_bonds + 1) - 1)
        ke[start + q] = 0.0
    end
    for jj in 1:n_bonds
        b = first_bond + jj - 1
        j = data.nlist[b]
        l2 = 0.0
        for a in 1:D
            l2 += (x[a, j] - x[a, i])^2
        end
        k = data.bond_damage[b] * data.volume[j] / (sqrt(l2) * l2)
        for cb in 1:D, ra in 1:D
            kab = k * (x[ra, j] - x[ra, i]) * (x[cb, j] - x[cb, i])
            ke[start + (cb - 1) * D + ra - 1] += kab
            ke[start + (jj * D + cb - 1) * D + ra - 1] -= kab
        end
    end
    return nothing
end

function neighbor_list(coordinates::AbstractMatrix{Float64}, horizon::Float64)
    n = size(coordinates, 2)
    return VectorOfVectors([[j
                             for j in 1:n
                             if j != i &&
                                norm(coordinates[:, j] - coordinates[:, i]) <=
                                horizon + 1e-12]
                            for i in 1:n])
end

"1D bar with `n` points, spacing 1, horizon 3."
function bar_1d(n::Int64; volume::Vector{Float64} = fill(1.0, n))
    coordinates = reshape(collect(0.0:(n - 1)), 1, n)
    nlist = neighbor_list(coordinates, 3.0)
    return (coordinates = coordinates, volume = volume,
            nlist_ptr = collect(Int64, internal_element_ptr(nlist)),
            nlist = collect(flatview(nlist)), bond_damage = ones(length(flatview(nlist))))
end

function bar_matrix(n::Int64)
    data = bar_1d(n; volume = [0.5 + 0.01i for i in 1:n])
    el = TestSpring{1}()
    pattern, K, ke = build_system(el, data, n)
    compute_elements!(ke, pattern, el, data)
    assemble!(K, pattern, ke)
    return K
end

"Direct solve with explicit sub-matrices; reference only."
function reference_solve(K::AbstractMatrix, prescribed::Vector{Int64},
                         u_p::Vector{Float64}, f::Vector{Float64})
    n = size(K, 1)
    free = setdiff(1:n, prescribed)
    u = zeros(n)
    u[prescribed] = u_p
    u[free] = K[free, free] \ (f[free] - K[free, prescribed] * u_p)
    return u
end
