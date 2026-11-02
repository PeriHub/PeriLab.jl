# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    run_perilab("bond_value_export", 1, true, folder_name)
    run_perilab("bond_value_export", 2, true, folder_name)
    run_perilab("bond_value_export_block_wise", 1, true, folder_name)
    run_perilab("bond_value_export_block_wise", 2, true, folder_name)
end
