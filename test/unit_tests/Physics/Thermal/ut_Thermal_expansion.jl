# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

include("../../../../src/Physics/Thermal/thermal_expansion.jl")
using .Thermal_expansion
#include("../../../../src/PeriLab.jl")
#import .PeriLab

@test Thermal_expansion.thermal_model_name() == "Thermal Expansion"

@testset "ut_thermal strain" begin
    dof = 2
    alpha = zeros(Int64, dof, dof)

    test_value = Thermal_expansion.thermal_strain(alpha, 0.0)
    @test test_value == zeros(dof, dof)
    alpha = zeros(Float64, dof, dof)

    test_value = Thermal_expansion.thermal_strain(alpha, 1)
    @test test_value == zeros(Float64, dof, dof)

    test_value = Thermal_expansion.thermal_strain(alpha, 0.0)
    @test test_value == zeros(Float64, dof, dof)
    alpha[1, 1] = -1
    alpha[2, 1] = 0
    alpha[1, 2] = 3
    alpha[2, 2] = 2
    temperature::Float64 = 1.5
    test_value = Thermal_expansion.thermal_strain(alpha, temperature)
    @test test_value[1:2, 1:2] == alpha .* temperature

end

@testset "ut_thermal_deformation" begin
    nnodes = 2
    dof = 2
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(2)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    temperature = test_data_manager.create_constant_node_field("Temperature", Float64, 1)
    undeformed_bond =
        test_data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    undeformed_bond_length =
        test_data_manager.create_constant_bond_field("Bond Length", Float64, 1)
    thermal_bond_deformation =
        test_data_manager.create_constant_bond_field("Thermal Deformation", Float64, dof)

    undeformed_bond[1][1, 1] = 0
    undeformed_bond[1][1, 2] = 1
    undeformed_bond_length[1][1] = 1
    undeformed_bond[1][2, 1] = 1
    undeformed_bond[1][2, 2] = 1
    undeformed_bond_length[1][2] = sqrt(2)

    undeformed_bond[2][1, 1] = -1
    undeformed_bond[2][1, 2] = -1
    undeformed_bond_length[2][1] = sqrt(2)
    undeformed_bond[2][2, 1] = 10
    undeformed_bond[2][2, 2] = -10
    undeformed_bond_length[2][2] = sqrt(200)
    undeformed_bond[2][3, 1] = 0.1
    undeformed_bond[2][3, 2] = 0
    undeformed_bond_length[2][3] = 0.1

    nodes = Vector{Int64}(1:nnodes)
    alpha = zeros(dof, dof)
    alpha[1, 1] = 1.0
    alpha[2, 2] = 1.0
    temperature .= 0
    thermal_bond_deformation = Thermal_expansion.thermal_deformation(
        nodes,
        alpha,
        temperature,
        undeformed_bond,
        thermal_bond_deformation,
    )
    for iID in nodes
        for jID in nn[iID]
            for i = 1:dof
                @test thermal_bond_deformation[iID][jID, i] == 0
            end
        end
    end
    temperature .= 1
    thermal_bond_deformation = Thermal_expansion.thermal_deformation(
        nodes,
        alpha,
        temperature,
        undeformed_bond,
        thermal_bond_deformation,
    )
    for iID in nodes
        for jID in nn[iID]
            @test thermal_bond_deformation[iID][jID, :] ==
                  .-undeformed_bond[iID][jID, 1:dof]
        end
    end

    temperature .= 2
    thermal_bond_deformation_test = Thermal_expansion.thermal_deformation(
        nodes,
        alpha,
        temperature,
        undeformed_bond,
        thermal_bond_deformation,
    )
    for iID in nodes
        for jID in nn[iID]
            @test isapprox(
                thermal_bond_deformation[iID][jID, :] .+ 1,
                .-undeformed_bond[iID][jID, 1:dof] .* 2 .+ 1,
            )
        end
    end

    temperature[1] = 2
    temperature[2] = -23

    alpha[1, 1] = -1.1
    alpha[2, 2] = 2.1
    thermal_bond_deformation = Thermal_expansion.thermal_deformation(
        nodes,
        alpha,
        temperature,
        undeformed_bond,
        thermal_bond_deformation,
    )
    for iID in nodes
        for jID in nn[iID]
            @test isapprox(
                thermal_bond_deformation[iID][jID, 1] + 1,
                .-undeformed_bond[iID][jID, 1] * alpha[1, 1] * temperature[iID] + 1,
            )
            @test isapprox(
                thermal_bond_deformation[iID][jID, 2] + 1,
                .-undeformed_bond[iID][jID, 2] * alpha[2, 2] * temperature[iID] + 1,
            )
        end
    end
end
