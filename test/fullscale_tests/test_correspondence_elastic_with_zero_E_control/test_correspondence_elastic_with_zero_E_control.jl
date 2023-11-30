# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../helper.jl")

folder_name = basename(@__FILE__)[1:end-3]
cd("fullscale_tests/" * folder_name) do
    # run_perilab("correspondence_strain_pstrain_zc_xx", 1, true, folder_name)
    run_perilab("correspondence_strain_pstrain_zc_xy", 1, true, folder_name)
    # run_perilab("correspondence_strain_pstrain_zc_yy", 1, true, folder_name)
    # run_perilab("correspondence_strain_pstress_zc_xx", 1, true, folder_name)
    run_perilab("correspondence_strain_pstress_zc_xy", 1, true, folder_name)
    # run_perilab("correspondence_strain_pstress_zc_yy", 1, true, folder_name)
end