# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../helper.jl")

folder_name = basename(@__FILE__)[1:end-3]
cd(folder_name) do
    run_perilab("displacements", 1, false, folder_name)
end