include("../../src/Support/data_manager.jl")
using Test
using MPI
using .Data_manager
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
testDatamanager.set_nnodes(nnodes)
nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
nn[:] = nneighbors
nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
@testset "init_distributed_lists rank $rank" begin
    for (id, list) in enumerate(nlist)
        @test length(list) == nneighbors[id]
    end
end

@testset "init_data rank $rank" begin

end
MPI.Finalize()