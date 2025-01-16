# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause



folder_name = basename(@__FILE__)[1:end-3]
cd("fullscale_tests/" * folder_name) do
    # run_perilab("strain_xx", 1, true, folder_name)
    run_perilab("strain_xy", 1, true, folder_name)
    # run_perilab("strain_xy", 3, true, folder_name)
    # run_perilab("strain_xz", 1, true, folder_name)
    # run_perilab("strain_yy", 1, true, folder_name)
    # run_perilab("strain_yz", 1, true, folder_name)
    # run_perilab("strain_zz", 1, true, folder_name)
end
