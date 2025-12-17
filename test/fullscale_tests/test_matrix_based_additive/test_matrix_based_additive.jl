# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    run_perilab("matrix_based_simple_printing", 1, true, folder_name)
    run_perilab("matrix_based_only_thermal", 1, true, folder_name)
    run_perilab("matrix_based_only_thermal_without_update", 1, true, folder_name)
end
