# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include(
    "../../../../../../src/Models/Material/Material_Models/BondBased/Bondbased_Elastic.jl",
)

using .Bondbased_Elastic
using Test
using TimerOutputs
#include("../../../../../../src/PeriLab.jl")
#using .PeriLab

#const to = TimerOutput()

@testset "material_name" begin
    @test Bondbased_Elastic.material_name() == "Bond-based Elastic"
    @test !(Bondbased_Elastic.fe_support())
end
@testset "compute_model" begin
    nodes = 2
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nodes)
    dof = 3
    test_data_manager.set_dof(dof)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    h = test_data_manager.create_constant_node_field("Horizon", Float64, 1)

    h[1:nodes] = 1:nodes
    bf = test_data_manager.create_constant_bond_field("Bond Forces", Float64, dof)

    nn[1] = 2
    nn[2] = 3

    bdN, bdNP1 = test_data_manager.create_bond_field("Bond Damage", Float64, 1)
    dbN, dbNP1 = test_data_manager.create_bond_field("Deformed Bond Geometry", Float64, dof)
    dbdN, dbdNP1 = test_data_manager.create_bond_field("Deformed Bond Length", Float64, 1)
    bg = test_data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    bd = test_data_manager.create_constant_bond_field("Bond Length", Float64, 1)
    for iID = 1:nodes
        bdNP1[iID] .= 1
        bd[iID] .= 1
        dbdNP1[iID] .= 1 + (-1)^iID * 0.1
        dbNP1[iID] .= 1
    end

    test_data_manager = Bondbased_Elastic.compute_model(
        test_data_manager,
        Vector{Int64}(1:nodes),
        Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0),
        1,
        0.0,
        0.0,
        to,
    )

    bf = test_data_manager.get_field("Bond Forces")
    @test isapprox(
        bf[1],
        [
            -0.13641852265019597 -0.13641852265019597 -0.13641852265019597
            -0.13641852265019597 -0.13641852265019597 -0.13641852265019597
        ],
    )
    @test isapprox(
        bf[2],
        [
            0.006975947180975936 0.006975947180975936 0.006975947180975936
            0.006975947180975936 0.006975947180975936 0.006975947180975936
            0.006975947180975936 0.006975947180975936 0.006975947180975936
        ],
    )


    bf[1][:, :] .= 0
    bf[2][:, :] .= 0

    test_data_manager = Bondbased_Elastic.compute_model(
        test_data_manager,
        Vector{Int64}(1:nodes),
        Dict(
            "Bulk Modulus" => 1.0,
            "Young's Modulus" => 1.0,
            "Symmetry" => "here is something",
        ),
        1,
        0.0,
        0.0,
        to,
    )

    bf = test_data_manager.get_field("Bond Forces")

    @test isapprox(
        bf[1],
        [
            -0.13641852265019597 -0.13641852265019597 -0.13641852265019597
            -0.13641852265019597 -0.13641852265019597 -0.13641852265019597
        ],
    )
    @test isapprox(
        bf[2],
        [
            0.006975947180975936 0.006975947180975936 0.006975947180975936
            0.006975947180975936 0.006975947180975936 0.006975947180975936
            0.006975947180975936 0.006975947180975936 0.006975947180975936
        ],
    )

    bf[1][:, :] .= 0
    bf[2][:, :] .= 0
    test_data_manager = Bondbased_Elastic.compute_model(
        test_data_manager,
        Vector{Int64}(1:nodes),
        Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0, "Symmetry" => "plane strain"),
        1,
        0.0,
        0.0,
        to,
    )

    bf = test_data_manager.get_field("Bond Forces")

    @test isapprox(
        bf[1],
        [
            -0.13058867125488846 -0.13058867125488846 -0.13058867125488846
            -0.13058867125488846 -0.13058867125488846 -0.13058867125488846
        ],
    )
    @test isapprox(
        bf[2],
        [
            0.01335565956015906 0.01335565956015906 0.01335565956015906
            0.01335565956015906 0.01335565956015906 0.01335565956015906
            0.01335565956015906 0.01335565956015906 0.01335565956015906
        ],
    )

    bf[1][:, :] .= 0
    bf[2][:, :] .= 0
    test_data_manager = Bondbased_Elastic.compute_model(
        test_data_manager,
        Vector{Int64}(1:nodes),
        Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0, "Symmetry" => "plane stress"),
        1,
        0.0,
        0.0,
        to,
    )

    bf = test_data_manager.get_field("Bond Forces")
    @test isapprox(bf[1][1, 1], -0.15915494309189532)
    @test isapprox(bf[1][1, 2], -0.15915494309189532)
    @test isapprox(bf[1][2, 1], -0.15915494309189532)
    @test isapprox(bf[1][2, 2], -0.15915494309189532)
    @test isapprox(bf[2][1, 1], 0.016277210088943856)
    @test isapprox(bf[2][1, 2], 0.016277210088943856)
    @test isapprox(bf[2][2, 1], 0.016277210088943856)
    @test isapprox(bf[2][2, 2], 0.016277210088943856)

end