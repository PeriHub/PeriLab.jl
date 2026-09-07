# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Reduction_template
using LinearAlgebra
using SparseArrays
export model_reduction_name
export reduce_matrices

function model_reduction_name()
    return "Reduction Template"
end

function reduce_matrices(K::AbstractMatrix{Float64}, M::Vector{Float64},
                         m::Vector{Int64}, s::Vector{Int64}, n_modes::Int64 = 1)
    return sparse(K), sparse(matrix(M))
end

end
