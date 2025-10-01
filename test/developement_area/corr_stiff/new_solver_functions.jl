# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
	compute_displacements!(K, non_BCs, u, F, F_temp, K_reduced, lu_fact, temp)

Compute displacements with prescribed displacement boundary conditions.
No memory allocations.

Arguments:
- K: Global stiffness matrix
- non_BCs: Free DOF indices
- u: Displacement vector (contains prescribed values at fixed DOFs)
- F: External force vector
- F_temp: Temporary force vector (pre-allocated)
- lu_fact: Pre-computed LU factorization of K[non_BCs, non_BCs]
- temp: Temporary vector for solution (length = number of free DOFs)
"""
function compute_displacements!(K::SparseMatrixCSC{Float64,Int64},
                                non_BCs::Vector{Int64},
                                u::Vector{Float64},
                                F::Vector{Float64},
                                F_temp::Vector{Float64},
                                temp::Vector{Float64})

    # Compute modified force: F_modified = F - K * u_prescribed
    # (u contains prescribed displacements at fixed DOFs)
    mul!(F_temp, K, u)  # F_temp = K * u (in-place, no allocation)
    F_temp .= F .- F_temp  # F_temp = F - K*u (in-place)

    # Solve for free DOFs: K_reduced * u_free = F_temp[non_BCs]
    #@views ldiv!(temp, lu(K[non_BCs, non_BCs]), F_temp[non_BCs])

    # Update displacement vector at free DOFs
    @views u[non_BCs] .= K[non_BCs, non_BCs] / F_temp[non_BCs]'

    return nothing
end
