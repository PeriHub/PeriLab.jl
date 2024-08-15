# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using PeriLab

folder_name = basename(@__FILE__)[1:end-3]
cd("fullscale_tests/" * folder_name) do
    PeriLab.main("dry_run.yaml"; silent = true, dry_run = true, output_dir = "temp")
    rm("temp", recursive = true)
end
