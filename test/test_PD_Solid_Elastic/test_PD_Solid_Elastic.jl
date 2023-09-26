# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../helper.jl")

cd(basename(@__FILE__)[1:end-3]) do
    run_perilab("strain_xx", 1, true)
    run_perilab("strain_xy", 1, true)
    run_perilab("strain_yy", 1, true)
end