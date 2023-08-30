include("../mesh_data.jl")
include("../../Support/data_manager.jl")
include("../../Support/Parameters/parameter_handling.jl")
import .Data_manager
using Test
import .Read_Mesh
@testset "ut_neighbors" begin
    nlist = []

    for i in 1:4
        append!(nlist, [collect(1:3*i*i-2)])
    end

    lenNlist = Read_Mesh.get_number_of_neighbornodes(nlist)

    for i in 1:4
        @test lenNlist[i] == 3 * i * i - 2
    end

end

@testset "ut_glob_to_loc" begin

    distribution = [1, 2, 3, 4, 5]
    glob_to_loc = Read_Mesh.glob_to_loc(distribution)
    len = length(distribution)
    #check trivial case of global and local are identical
    for id in 1:len
        @test distribution[id] == glob_to_loc[id]
    end
    @test length(distribution) == length(glob_to_loc)
    # reverse -> glob_to_loc_to_glob
    distribution = [1, 4, 2, 5, 6]
    glob_to_loc = Read_Mesh.glob_to_loc(distribution)
    for id in 1:len
        @test distribution[id] == distribution[glob_to_loc[distribution[id]]]
    end

end

@testset "ut_define_nsets" begin

    numbers = [11, 12, 13, 44, 125]
    lenNumbers = length(numbers)
    filename = "test.txt"
    file = open(filename, "w")
    for number in numbers
        println(file, number)
    end
    close(file)
    params = Dict("Discretization" => Dict("Node Sets" => Dict("Nset_1" => "1 2 3 4 5 6 7", "Nset_2" => filename)))
    testDatamanager = Data_manager
    @test testDatamanager.get_nnsets() == 0
    Read_Mesh.define_nsets(params, testDatamanager)
    @test testDatamanager.get_nnsets() == 2
    nsets = testDatamanager.get_nsets()
    @test nsets["Nset_1"] == [1, 2, 3, 4, 5, 6, 7]
    @test nsets["Nset_2"] == [11, 12, 13, 44, 125]

    rm(filename)
end