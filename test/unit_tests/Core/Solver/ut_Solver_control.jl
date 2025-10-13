# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using MPI
using TimerOutputs

@testset "ut_get_block_nodes" begin
    block_ids = [1, 1, 1, 2, 2, 3, 3, 3, 3, 1, 1, 2, 3, 3, 1, 1, 2]
    block_nodes = PeriLab.Solver_Manager.get_block_nodes(block_ids, length(block_ids))
    @test block_nodes[1] == [1, 2, 3, 10, 11, 15, 16]
    @test block_nodes[2] == [4, 5, 12, 17]
    @test block_nodes[3] == [6, 7, 8, 9, 13, 14]
end
