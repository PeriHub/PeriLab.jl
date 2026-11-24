# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# 2D bar, block 2 (x in 2..8 mm) condensed, displacement load on the column x = 4.8 mm in
# the middle of it, kept by "Retained Node Sets". The reduced runs are compared with the
# unreduced one on the nodes of block 1 at the end time.

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    function last_displacements(name)
        exo = ExodusDatabase(name * ".e", "r")
        n = read_number_of_time_steps(exo)
        u = hcat(read_values(exo, NodalVariable, n, 1, "Displacementsx"),
                 read_values(exo, NodalVariable, n, 1, "Displacementsy"))
        close(exo)
        rm(name * ".e")
        return u
    end
    block_1 = [i for i in 1:156 if !(0.0019 < ((i - 1) % 26) * 0.0004 < 0.0081)]
    PeriLab.run("full.yaml"; silent = true)
    full = last_displacements("full")[block_1, :]
    # Craig-Bampton with all modes and without material points is exact up to the
    # prescribed acceleration taken one step late; material points add the mismatch
    # between their forces and the stiffness matrix at the interface
    for (name, tolerance) in [("cb_nomp", 1e-3), ("cbc_nomp", 1e-3), ("cb", 3e-2),
            ("sc", 3e-2)]
        @test PeriLab.run(name * ".yaml"; silent = true) == 0
        @test maximum(abs, last_displacements(name)[block_1, :] - full) <=
              tolerance * maximum(abs, full)
    end
    # the load column is condensed without the option: rejected
    @test PeriLab.run("no_option.yaml"; silent = true) == 1
end
