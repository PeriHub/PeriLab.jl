# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../helper.jl")

folder_name = dirname(@__FILE__)
run_mpi_test("ut_MPI.jl", 3, true, folder_name)
