# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    # a checkpoint left over from a previous test run would make reload1 start at its end time
    rm("restart"; recursive = true, force = true)
    run_perilab("reload1", 1, false, folder_name; reload = true)
    run_perilab("reload2", 1, true, folder_name; reload = true)
end
