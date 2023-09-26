# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../src/Support/data_manager.jl")
include("../../src/IO/mesh_data.jl")
using Test
using MPI
using .Data_manager
using .Read_Mesh
include("../../src/MPI_communication/MPI_communication.jl")
MPI.Init()
testDatamanager = Data_manager
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    nnodes = 3
    nneighbors = [2, 3, 2]
else
    nnodes = 4
    nneighbors = [2, 3, 2, 5]
end
testDatamanager.set_nmasters(nnodes)
nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
nn[:] = nneighbors
nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
@testset "init_distributed_lists rank $rank" begin
    for (id, list) in enumerate(nlist)
        @test length(list) == nneighbors[id]
    end
end
distribution = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [7, 8, 9, 10, 11, 12, 3, 4, 5, 6]]

if rank == 0
    send_msg = [2, 3, 3, 2, 3, 4, 4, 3, 2, 3, 3, 2]
end
glob_to_loc = Read_Mesh.glob_to_loc(distribution)
@testset "init_data rank $rank" begin
    lenNlist[:] = send_vector_from_root_to_core_i(comm, send_msg, lenNlist, distribution)
    nlistCore = datamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)

end
MPI.Finalize()