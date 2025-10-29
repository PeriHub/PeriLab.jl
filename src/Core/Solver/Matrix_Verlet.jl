# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# craig_bambton

function craig_bampton(K, M, m, s, n_modes)
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
