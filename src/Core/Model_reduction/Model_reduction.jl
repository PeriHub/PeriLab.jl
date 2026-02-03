# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Model_reduction
using LinearAlgebra
using SparseArrays
function guyan_reduction(K::AbstractMatrix{Float64}, M_diag::Vector{Float64},
                         m::Vector{Int64}, s::Vector{Int64})
    nm = length(m)
    ns = length(s)

    # Extract submatrices
    K_mm = K[m, m]
    K_ms = K[m, s]
    K_ss = K[s, s]
    K_sm = K[s, m]

    # Compute transformation matrix: T = -K_ss^(-1) * K_sm
    K_ss_fact = lu(K_ss)
    K_sm_dense = Matrix(K_sm)
    T = similar(K_sm_dense)
    ldiv!(T, K_ss_fact, K_sm_dense)
    T .*= -1.0

    # Stiffness reduction: K_reduced = K_mm + K_ms * T
    temp_K = zeros(nm, nm)
    mul!(temp_K, K_ms, T)
    K_reduced = K_mm + temp_K

    # Mass reduction: M_reduced = M_mm + T^T * M_ss * T
    # Step 1: temp_M = Diagonal(M_ss) * T
    M_ss_diag = Diagonal(M_diag[s])
    temp_M = zeros(ns, nm)
    mul!(temp_M, M_ss_diag, T)  # Diagonal * Matrix multiplication

    # Step 2: M_reduced = T^T * temp_M
    M_reduced = zeros(nm, nm)
    mul!(M_reduced, T', temp_M)

    # Step 3: Add M_mm diagonal
    @inbounds for i in 1:nm
        M_reduced[i, i] += M_diag[m[i]]
    end

    return sparse(K_reduced), sparse(M_reduced)
end

#function stiffness_reduction(K::AbstractMatrix{Float64}, m::Vector{Int64}, s::Vector{Int64})
#    return K[m, m] - K[m, s] / Matrix(K[s, s]) * K[s, m]
#end
#
#function reduce_mass_consistent(M_diag, K, master_nodes, slave_nodes)
#    # Transformationsmatrix: u_slave = T * u_master
#    # T = -K_ss^{-1} * K_sm
#    Ksm = K[slave_nodes, master_nodes]
#    Kss = K[slave_nodes, slave_nodes]
#    T = -Kss \ Matrix(Ksm)
#
#    # Massenreduktion: M_red = M_mm + T^T * M_ss * T
#    M_mm = Diagonal(M_diag[master_nodes])
#    M_ss = Diagonal(M_diag[slave_nodes])
#
#    M_reduced = M_mm + T' * M_ss * T
#
#    return sparse(M_mm)
#end

function craig_bampton(K::AbstractMatrix{Float64}, M::AbstractMatrix{Float64},
                       m::Vector{Int64}, s::Vector{Int64}, n_modes::Int64)
    # Partitionierung
    K_mm = K[m, m]
    K_ms = K[m, s]
    K_ss = K[s, s]
    K_sm = K[s, m]

    M_mm = M[m, m]
    M_ms = M[m, s]
    M_ss = M[s, s]
    M_sm = M[s, m]

    # 1. CONSTRAINT MODES (statisch, wie bei Guyan)
    Φ_c = -(K_ss \ K_sm)  # Constraint modes

    # 2. FIXED-INTERFACE NORMAL MODES (neu!)
    # Eigenwertproblem mit fixierten Master-DOFs
    eigenvals, eigenvecs = eigen(K_ss, M_ss)

    # Nur die ersten n_modes behalten
    Φ_n = eigenvecs[:, 1:n_modes]  # Normal modes
    λ_n = eigenvals[1:n_modes]      # Eigenvalues

    # 3. TRANSFORMATIONSMATRIX
    # [x_m]   [I    0  ] [x_m  ]
    # [x_s] = [Φ_c  Φ_n] [η_n  ]
    #
    # x_m: physikalische Master-DOFs
    # η_n: modale Koordinaten

    T_full = [Matrix(I, length(m), length(m)) zeros(length(m), n_modes);
              Φ_c Φ_n]

    # 4. REDUZIERTE MATRIZEN
    K_craig = T_full' * [K_mm K_ms; K_sm K_ss] * T_full
    M_craig = T_full' * [M_mm M_ms; M_sm M_ss] * T_full

    # Die reduzierten Matrizen haben spezielle Struktur:
    # K_craig = [K_mm + K_ms*Φ_c    K_ms*Φ_n  ]
    #           [Φ_n'*K_sm          Λ_n       ]
    #
    # M_craig = [M_mm + M_ms*Φ_c + ...    ...  ]
    #           [...                     I     ]

    return K_craig, M_craig, Φ_c, Φ_n, λ_n
end

end
