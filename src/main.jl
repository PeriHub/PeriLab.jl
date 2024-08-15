# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Pkg
using MPI
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
    Pkg.instantiate()
end

MPI.Init()
import PeriLab
PeriLab.main()
