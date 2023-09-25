include("../../../../../../src/Physics/Material/Material_Models/Ordinary/Ordinary.jl")
using Test
using .Ordinary
@testset "compute_weighted_volume" begin
    weightedTest = zeros(9)
    # from Peridigm
    weightedTest[1] = 6.0311182727183619
    weightedTest[2] = 15.723986925301443
    weightedTest[3] = 14.862398600627394
    weightedTest[4] = 14.86239860062739
    weightedTest[5] = 5.3849270292128226
    weightedTest[6] = 7.9696920032349787
    weightedTest[7] = 12.277633626605239
    weightedTest[8] = 6.2465153538868741
    weightedTest[9] = 19.601134386334675

    # values taken from the variable manager in Debug
    # comparison between Peridigm and PeriLab
    nnodes = 9
    nneighbors = [3, 5, 6, 4, 4, 5, 6, 4, 7]

    nlist = Any[Any[2, 3, 4], Any[1, 3, 4, 7, 9], Any[1, 2, 4, 6, 7, 9], Any[1, 2, 3, 9], Any[6, 7, 8, 9], Any[3, 5, 7, 8, 9], Any[2, 3, 5, 6, 8, 9], Any[5, 6, 7, 9], Any[2, 3, 4, 5, 6, 7, 8]]

    bond_geometry = Any[Float32[0.0 1.0 1.0; 1.0 1.0 1.4142135; 2.0 0.0 2.0], Float32[0.0 -1.0 1.0; 1.0 0.0 1.0; 2.0 -1.0 2.236068; 1.0 2.0 2.236068; 2.0 1.5 2.5], Float32[-1.0 -1.0 1.4142135; -1.0 0.0 1.0; 1.0 -1.0 1.4142135; 1.0 2.0 2.236068; 0.0 2.0 2.0; 1.0 1.5 1.8027756], Float32[-2.0 0.0 2.0; -2.0 1.0 2.236068; -1.0 1.0 1.4142135; 0.0 2.5 2.5], Float32[0.0 -1.0 1.0; -1.0 -1.0 1.4142135; -1.0 0.0 1.0; 0.0 -1.5 1.5], Float32[-1.0 -2.0 2.236068; 0.0 1.0 1.0; -1.0 0.0 1.0; -1.0 1.0 1.4142135; 0.0 -0.5 0.5], Float32[-1.0 -2.0 2.236068; 0.0 -2.0 2.0; 1.0 1.0 1.4142135; 1.0 0.0 1.0; 0.0 1.0 1.0; 1.0 -0.5 1.118034], Float32[1.0 0.0 1.0; 1.0 -1.0 1.4142135; 0.0 -1.0 1.0; 1.0 -1.5 1.8027756], Float32[-2.0 -1.5 2.5; -1.0 -1.5 1.8027756; 0.0 -2.5 2.5; 0.0 1.5 1.5; 0.0 0.5 0.5; -1.0 0.5 1.118034; -1.0 1.5 1.8027756]]

    bond_damage = Any[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]

    omega = Float32[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    volume = Float32[0.8615883, 0.8615883, 0.8615883, 0.8615883, 0.8615883, 0.8615883, 0.8615883, 0.8615883, 0.8615883]

    weighted_volume = Ordinary.compute_weighted_volume(1:nnodes, nneighbors, nlist, bond_geometry, bond_damage, omega, volume)

    for iID in 1:nnodes
        @test weighted_volume[iID] / weightedTest[iID] - 1 < 1e-6
    end

end

nnodes = 2
nneighbors = [1, 1]

bond_geometry = [Matrix{Float32}(undef, 1, 3), Matrix{Float32}(undef, 1, 3)]
bond_geometry[1][1] = 1
bond_geometry[1][2] = 0
bond_geometry[1][3] = 1
bond_geometry[2][1] = -1
bond_geometry[2][2] = 0
bond_geometry[2][3] = 1

deformed_bond = [Matrix{Float32}(undef, 1, 3), Matrix{Float32}(undef, 1, 3)]
deformed_bond[1][1, 1] = 2
deformed_bond[1][1, 2] = 0
deformed_bond[1][1, 3] = 2
deformed_bond[2][1, 1] = -2
deformed_bond[2][1, 2] = 0
deformed_bond[2][1, 3] = 2
bond_damage = [Vector{Float32}(undef, 1), Vector{Float32}(undef, 1)]
bond_damage[1][1] = 1
bond_damage[2][1] = 1

volume = ones(Float32, 2)
weighted_volume = ones(Float32, 2)
omega = ones(Float32, 2)

nlist = [Vector{Int64}(undef, 1), Vector{Int64}(undef, 1)]
nlist[1][1] = 2
nlist[2][1] = 2
@testset "compute_dilatation" begin
    theta = Ordinary.compute_dilatation(1:nnodes, nneighbors, nlist, bond_geometry, deformed_bond, bond_damage, volume, weighted_volume, omega)
    @test theta[1] == 3.0
    @test theta[2] == 3.0
end