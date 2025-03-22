# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../../src/Models/Material/Material_Models/Ordinary/Ordinary.jl")
using Test
using .Ordinary
@testset "compute_weighted_volume" begin
    weighted_volume = zeros(9)
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

    nlist = [
        [2, 3, 4],
        [1, 3, 4, 7, 9],
        [1, 2, 4, 6, 7, 9],
        [1, 2, 3, 9],
        [6, 7, 8, 9],
        [3, 5, 7, 8, 9],
        [2, 3, 5, 6, 8, 9],
        [5, 6, 7, 9],
        [2, 3, 4, 5, 6, 7, 8]
    ]

    undeformed_bond_length = [
        Float64[1.0; 1.4142135; 2.0],
        Float64[1.0; 1.0; 2.236068; 2.236068; 2.5],
        Float64[1.4142135; 1.0; 1.4142135; 2.236068; 2.0; 1.8027756],
        Float64[2.0; 2.236068; 1.4142135; 2.5],
        Float64[1.0; 1.4142135; -1.0; 1.5],
        Float64[2.236068; 1.0; -1.0; 1.4142135; 0.5],
        Float64[2.236068; 2.0; 1.4142135; 1.0; 1.0; 1.118034],
        Float64[1.0; 1.4142135; 1.0; 1.8027756],
        Float64[2.5; 1.8027756; 2.5; 1.5; 0.5; 1.118034; 1.8027756]
    ]

    bond_damage = [
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    ]

    omega = [
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    ]

    volume = Float64[0.8615883,
                     0.8615883,
                     0.8615883,
                     0.8615883,
                     0.8615883,
                     0.8615883,
                     0.8615883,
                     0.8615883,
                     0.8615883]
    vec = Vector{Int64}(1:nnodes)
    Ordinary.compute_weighted_volume(weighted_volume,
                                     vec,
                                     nlist,
                                     undeformed_bond_length,
                                     bond_damage,
                                     omega,
                                     volume)

    for iID in 1:nnodes
        @test weighted_volume[iID] / weightedTest[iID] - 1 < 1e-6
    end
end

nnodes = 2
nneighbors = [1, 1]

undeformed_bond_length = [Vector{Float64}(undef, 1), Vector{Float64}(undef, 1)]
undeformed_bond_length[1][1] = 1
undeformed_bond_length[2][1] = 1

deformed_bond_length = [Vector{Float64}(undef, 1), Vector{Float64}(undef, 1)]
deformed_bond_length[1][1] = 2
deformed_bond_length[2][1] = 2
bond_damage = [Vector{Float64}(undef, 1), Vector{Float64}(undef, 1)]
bond_damage[1][1] = 1
bond_damage[2][1] = 1

volume = ones(Float64, 2)
weighted_volume = ones(Float64, 2)
theta = zeros(Float64, 2)
omega = [ones(Float64, 2), ones(Float64, 2)]

nlist = [Vector{Int64}(undef, 1), Vector{Int64}(undef, 1)]
nlist[1][1] = 2
nlist[2][1] = 2
@testset "compute_dilatation" begin
    vec = Vector{Int64}(1:nnodes)
    theta = zeros(Float64, 2)
    Ordinary.compute_dilatation(vec,
                                nneighbors,
                                nlist,
                                undeformed_bond_length,
                                deformed_bond_length,
                                bond_damage,
                                volume,
                                weighted_volume,
                                omega,
                                theta)
    @test theta[1] == 3.0
    @test theta[2] == 3.0
    weighted_volume[1] = 0
    Ordinary.compute_dilatation(vec,
                                nneighbors,
                                nlist,
                                undeformed_bond_length,
                                deformed_bond_length,
                                bond_damage,
                                volume,
                                weighted_volume,
                                omega,
                                theta)
    @test theta[1] == 0.0
    @test theta[2] == 3.0
end

@testset "calculate_symmetry_params" begin
    @test Ordinary.calculate_symmetry_params("3D", 1.0, 1.0) == (15, 1, 3)
    test_val = Ordinary.calculate_symmetry_params("plane stress", 1.0, 1.0)
    @test test_val[1] == 8
    @test isapprox(test_val[2], 0.5714285714285714)
    @test isapprox(test_val[3], 0.5714285714285714)
    @test Ordinary.calculate_symmetry_params("plane strain", 1.0, 1.0) == (8, 2 / 3, 8 / 9)
end

@testset "ut_get_bond_forces" begin
    vec = Vector{Int64}(1:nnodes)
    dof = 2
    nBonds = fill(dof, nnodes)
    bond_force_length = [fill(1.0, n) for n in nBonds]
    deformed_bond = [[fill(0.0, dof) for j in 1:n] for n in nBonds]
    deformed_bond[1][1][1] = 1
    deformed_bond[1][2][1] = 1
    bond_force = [[fill(0.0, dof) for j in 1:n] for n in nBonds]
    temp = [fill(0.0, n) for n in nBonds]
    bond_force = Ordinary.get_bond_forces(vec,
                                          bond_force_length,
                                          deformed_bond,
                                          deformed_bond_length,
                                          bond_force,
                                          temp)
    @test bond_force == [[[0.5, 0.0], [0.5, 0.0]], [[0.0, 0.0], [0.0, 0.0]]]
    # deformed_bond_length[2][1] = 0
    # @test isnothing(
    #     Ordinary.get_bond_forces(
    #         vec,
    #         bond_force_length,
    #         deformed_bond,
    #         deformed_bond_length,
    #         bond_force,
    #         temp,
    #     ),
    # )
end
