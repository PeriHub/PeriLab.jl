# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    run_perilab("penalty_contact", 1, true, folder_name)
    run_perilab("penalty_contact_options", 1, true, folder_name)
    run_perilab("penalty_contact_friction", 1, true, folder_name)
    run_perilab("penalty_contact_options", 2, true, folder_name)
    # run_perilab("penalty_contact_options", 3, true, folder_name)
end
