# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

# Run with 3 ranks: mpiexec -n 3 julia --project test/multifunctional_test/test_MPI_matrix.jl
using Test
using MPI
using Logging
import PeriLab # not `using`: it would run PeriLab.main at exit
using PeriLab.Solver_Manager.PD_Matrix: create_halo_exchange, exchange!,
                                        accumulate_ghosts!
using PeriLab.Solver_Manager.Matrix_Adapter: perilab_halo
Logging.disable_logging(Logging.Warn)

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
n_ranks = MPI.Comm_size(comm)

@testset "rank $rank" begin
    if n_ranks == 3
        @testset "ghosts in arbitrary owner order" begin
            dof = 2
            own = collect((4rank + 1):(4rank + 4))
            # descending owners and interleaved points: not grouped like the receive order
            ghosts = [4r + k for r in 2:-1:0 if r != rank for k in (3, 1)]
            owner = [(g - 1) ÷ 4 for g in ghosts]
            @test !issorted(owner)
            halo = create_halo_exchange(comm, own, ghosts, owner, dof)
            x = zeros(dof * (length(own) + length(ghosts)))
            for (i, g) in enumerate(own), c in 1:dof
                x[(i - 1) * dof + c] = 10g + c
            end
            exchange!(halo, x)
            n_own = length(own)
            @test all(x[(n_own + i - 1) * dof + c] == 10g + c
                      for (i, g) in enumerate(ghosts) for c in 1:dof)
        end

        @testset "three ranks, end ranks" begin
            # 1D chain, rank r owns points 2r+1, 2r+2; ghosts are the direct neighbours
            own = [2rank + 1, 2rank + 2]
            ghosts = rank == 0 ? [3] : rank == 1 ? [5, 2] : [4]
            owner = [(g - 1) ÷ 2 for g in ghosts]
            halo = create_halo_exchange(comm, own, ghosts, owner, 1)
            x = vcat(Float64.(own), zeros(length(ghosts)))
            exchange!(halo, x)
            @test x[3:end] == Float64.(ghosts)
            y = vcat(zeros(2), ones(length(ghosts)))
            accumulate_ghosts!(halo, y)
            # point 2 (rank 0) and 5 (rank 2) are ghosts of rank 1, points 3 and 4 (rank 1)
            # are ghosts of rank 0 and rank 2
            expected = rank == 0 ? [0.0, 1.0] : rank == 1 ? [1.0, 1.0] : [1.0, 0.0]
            @test y[1:2] == expected
        end
    end

    @testset "perilab_halo" begin
        folder = joinpath(@__DIR__, "..", "fullscale_tests", "test_DCB")
        PeriLab.Data_Manager.initialize_data()
        PeriLab.IO.initialize_data(joinpath(folder, "DCBmodel_correspondence.yaml"), folder,
                                   comm)
        DM = PeriLab.Data_Manager
        nnodes = DM.get_nnodes()
        n_ghost = DM.get_num_responder()
        @test n_ghost > 0
        coor = DM.get_field("Coordinates")
        dof = DM.get_dof()
        halo = perilab_halo(dof)
        x = zeros(dof * (nnodes + n_ghost))
        for i in 1:nnodes, c in 1:dof
            x[(i - 1) * dof + c] = coor[i, c]
        end
        exchange!(halo, x)
        @test all(x[(i - 1) * dof + c] == coor[i, c]
                  for i in (nnodes + 1):(nnodes + n_ghost) for c in 1:dof)
    end
end
MPI.Finalize()
