using Test
include("../Solver_control.jl")
import .Solver
@testset "ut_get_blockNodes" begin
    blockIDs = [1, 1, 1, 2, 2, 3, 3, 3, 3, 1, 1, 2, 3, 3, 1, 1, 2]
    blockNodes = Solver.get_blockNodes(blockIDs)
    @test blockNodes[1] == [1, 2, 3, 10, 11, 15, 16]
    @test blockNodes[2] == [4, 5, 12, 17]
    @test blockNodes[3] == [6, 7, 8, 9, 13, 14]
end