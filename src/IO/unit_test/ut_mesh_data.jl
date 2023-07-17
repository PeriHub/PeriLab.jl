include("../mesh_data.jl")
using Test

@testset "ut_neighbors" begin
    nlist = []

    for i in 1:4
        append!(nlist, [collect(1:3*i*i-2)])
    end

    lenNlist = get_number_of_neighbornodes(nlist)

    for i in 1:4
        @test lenNlist[i] == 3 * i * i - 2
    end




end