# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    run_perilab("linear_static_matrix_based_solver_2D", 1, true, folder_name)
    run_perilab("linear_static_matrix_based_solver_3D", 1, true, folder_name)
    # the DCB decks have no x constraint: the x displacement contains a rigid body
    # translation that only depends on round-off, so it is not compared
    run_perilab("linear_static_matrix_based_solver_DCBmodel", 1, true, folder_name * "_DCB")
    run_perilab("linear_static_matrix_based_solver_switch_DCBmodel", 1, true,
                folder_name * "_DCB")
    run_perilab("linear_static_matrix_based_solver_fallback_DCBmodel", 1, true,
                folder_name * "_DCB")
    run_perilab("linear_static_matrix_based_tension", 1, true, folder_name)
    run_perilab("linear_static_matrix_based_thermal_expansion", 1, true, folder_name)
    run_perilab("linear_static_matrix_based_bending", 1, true, folder_name)
    run_perilab("linear_static_matrix_based_shear", 1, true, folder_name)
end
