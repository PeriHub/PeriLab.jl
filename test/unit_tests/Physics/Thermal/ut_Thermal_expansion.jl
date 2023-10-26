# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

include("../../../../src/Physics/Thermal/Thermal_expansion.jl")
using .Thermal_expansion
include("../../../../src/Support/data_manager.jl")

@testset "ut_thermal_deformation" begin
    nnodes = 2
    dof = 2
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(2)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    temperature = testDatamanager.create_constant_node_field("Temperature", Float64, 1)
    undeformed_bond = testDatamanager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    thermal_bond_deformation = testDatamanager.create_constant_bond_field("Thermal Deformation", Float64, dof + 1)

    undeformed_bond[1][1, 1] = 0
    undeformed_bond[1][1, 2] = 1
    undeformed_bond[1][1, 3] = 1
    undeformed_bond[1][1, 1] = 1
    undeformed_bond[1][1, 2] = 1
    undeformed_bond[1][1, 3] = sqrt(2)

    undeformed_bond[2][1, 1] = -1
    undeformed_bond[2][1, 2] = -1
    undeformed_bond[2][1, 3] = sqrt(2)
    undeformed_bond[2][2, 1] = 10
    undeformed_bond[2][2, 2] = -10
    undeformed_bond[2][2, 3] = sqrt(200)
    undeformed_bond[2][3, 1] = 0.1
    undeformed_bond[2][3, 2] = 0
    undeformed_bond[2][3, 3] = 0.1

    nodes = Vector{Int64}(1:nnodes)
    alpha = 1.0
    temperature .= 0
    thermal_bond_deformation = Thermal_expansion.thermal_deformation(nodes, alpha, temperature, undeformed_bond, thermal_bond_deformation)
    for iID in nodes
        for jID in nn[iID]
            for i in 1:dof+1
                @test thermal_bond_deformation[iID][jID, i] == 0
            end
        end
    end

    temperature .= 1
    thermal_bond_deformation_test = Thermal_expansion.thermal_deformation(nodes, alpha, temperature, undeformed_bond, thermal_bond_deformation)
    for iID in nodes
        for jID in nn[iID]
            @test isapprox(thermal_bond_deformation_test[iID][jID, :] .+ 1, undeformed_bond[iID][jID, :] .+ 1)
        end
    end
    temperature .= 2
    thermal_bond_deformation_test = Thermal_expansion.thermal_deformation(nodes, alpha, temperature, undeformed_bond, thermal_bond_deformation)
    for iID in nodes
        for jID in nn[iID]
            @test isapprox(thermal_bond_deformation_test[iID][jID, :] .+ 1, undeformed_bond[iID][jID, :] .* 2 .+ 1)
        end
    end
    temperature[1] = 2
    temperature[1] = -23
    alpha = -1.1
    thermal_bond_deformation_test = Thermal_expansion.thermal_deformation(nodes, alpha, temperature, undeformed_bond, thermal_bond_deformation)
    for iID in nodes
        for jID in nn[iID]
            @test isapprox(thermal_bond_deformation_test[iID][jID, :] .+ 1, undeformed_bond[iID][jID, :] .* (-1.1) * temperature[iID] .+ 1)
        end
    end
end
@testset "ut_thermal strain" begin
    nnodes = 2
    dof = 2
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(2)
    temperature = testDatamanager.create_constant_node_field("Temperature", Float64, 1)
    thermal_strain_tensor = testDatamanager.create_constant_node_field("Thermal Strain", Float64, "Matrix", dof)
    nodes = Vector{Int64}(1:nnodes)

    alpha = zeros(dof, dof)
    temperature .= 0.0

    test_value = Thermal_expansion.thermal_strain(nodes, alpha, temperature, thermal_strain_tensor)
    @test test_value == thermal_strain_tensor
    temperature .= 1.0
    test_value = Thermal_expansion.thermal_strain(nodes, alpha, temperature, thermal_strain_tensor)
    @test test_value == thermal_strain_tensor
    alpha = ones(dof, dof)
    temperature .= 0.0
    test_value = Thermal_expansion.thermal_strain(nodes, alpha, temperature, thermal_strain_tensor)
    @test test_value == thermal_strain_tensor
    alpha[1, 1] = -1
    alpha[2, 1] = 0
    alpha[1, 2] = 3
    alpha[2, 2] = 2
    temperature[1] = 1.5
    temperature[2] = 1.0
    test_value = Thermal_expansion.thermal_strain(nodes, alpha, temperature, thermal_strain_tensor)
    for iID in nodes
        @test test_value[iID, 1:2, 1:2] == alpha .* temperature[iID]
    end
    @test test_value == thermal_strain_tensor


end

