# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../helper.jl")

cd(basename(@__FILE__)[1:end-3]) do
    run_perilab("correspondence_strain_pstrain_xx_hg", 1, true)
    run_perilab("correspondence_strain_pstrain_xy_hg", 1, true)
    run_perilab("correspondence_strain_pstrain_yy_hg", 1, true)
    run_perilab("correspondence_strain_pstress_xx_hg", 1, true)
    run_perilab("correspondence_strain_pstress_xy_hg", 1, true)
    run_perilab("correspondence_strain_pstress_yy_hg", 1, true)
end