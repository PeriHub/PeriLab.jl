include("../exodus_export.jl")
include("../../Support/data_manager.jl")
using Test
import .Data_manager
import .Write_Exodus_Results
if !isdir("tmp")
    mkdir("tmp")
end
filename = "./tmp/" * "test.e"
@testset "ut_create_result_file" begin
    nnodes = 4
    dof = 3
    exo = Write_Exodus_Results.create_result_file(filename, nnodes, dof, 1, 0)
    @test isfile(filename)
    close(exo)
    rm(filename)
end
