# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../src/Core/Data_manager.jl")
include("../../src/IO/mesh_data.jl")
using Test
using MPI
using .Data_Manager
using .MPI_Communication: send_vector_from_root_to_core_i

include("../../src/MPI_communication/MPI_communication.jl")
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    nnodes = 3
    nneighbors = [2, 3, 2]
else
    nnodes = 4
    nneighbors = [2, 3, 2, 5]
end
Data_Manager.set_num_controller(nnodes)
nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
nn = nneighbors
nlist = Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
@testset "init_distributed_lists rank $rank" begin
    for (id, list) in enumerate(nlist)
        @test length(list) == nneighbors[id]
    end
end
distribution = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [7, 8, 9, 10, 11, 12, 3, 4, 5, 6]]

if rank == 0
    send_msg = [2, 3, 3, 2, 3, 4, 4, 3, 2, 3, 3, 2]
end
create_global_to_local_mapping = create_global_to_local_mapping(distribution)
@testset "init_data rank $rank" begin
    length_nlist = send_vector_from_root_to_core_i(comm, send_msg, length_nlist,
                                                   distribution)
    nlist_core = Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
end
MPI.Finalize()
