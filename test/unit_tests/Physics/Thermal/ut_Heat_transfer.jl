# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

include("../../../../src/Physics/Thermal/heat_transfer.jl")
using .Heat_transfer
# include("../../../../src/Support/data_manager.jl")

@test Heat_transfer.thermal_model_name() == "Heat Transfer"

@testset "ut_calculate_specific_volume" begin
    nnodes = 10
    dof = 2
    nodes = Vector{Int64}(1:nnodes)
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.clear_data_manager()
    test_Data_manager.set_num_controller(nnodes)
    test_Data_manager.set_dof(dof)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 3
    nn[5] = 4
    nn[6] = 3
    nn[7] = 2
    nn[8] = 3
    nn[9] = 2
    nn[10] = 1
    nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [4, 2]
    nlist[2] = [1, 3, 5]
    nlist[3] = [2, 6]
    nlist[4] = [1, 5, 7]
    nlist[5] = [2, 4, 6, 8]
    nlist[6] = [3, 5, 9]
    nlist[7] = [4, 8]
    nlist[8] = [5, 7, 9]
    nlist[9] = [6, 8]
    nlist[10] = [9]
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
    volume .= 0.25
    specific_volume = test_Data_manager.create_constant_node_field("specific_volume", Float64, 1)
    horizon = test_Data_manager.create_constant_node_field("Horizon", Float64, 1)
    horizon .= 0.55
    active = fill(true, nnodes)
    result = Heat_transfer.calculate_specific_volume(nodes, nlist, volume, active, specific_volume, dof, horizon)
    @test result == [1.0, 0.6666666666666666, 1.0, 0.6666666666666666, 0.5, 0.6666666666666666, 1.0, 0.6666666666666666, 1.0, 2.0]
end