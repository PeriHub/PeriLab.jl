# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    run_perilab("critical_energy_tension", 1, true, folder_name)
    run_perilab("critical_energy_interface", 1, true, folder_name)
    run_perilab("critical_energy_pressure", 1, true, folder_name)
    run_perilab("critical_energy_field", 1, true, folder_name)
end
