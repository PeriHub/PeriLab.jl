# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../src/Physics/Material/Material_Models/Correspondence.jl")
#include("../../../../../src/PeriLab.jl")
# include("../../../../../src/Core/data_manager.jl")
using Test
using .Correspondence
#using .PeriLab
@testset "zero_energy_mode_compensation_exception" begin

    test_data_manager = PeriLab.Data_manager
    nnodes = 2
    test_data_manager.set_num_controller(nnodes)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nodes = Vector{Int64}(1:2)

    @test Correspondence.zero_energy_mode_compensation(test_data_manager, nodes, Dict(), 0.0, 0.0) == test_data_manager
end
@testset "rotate_second_order_tensor" begin

    angles = [0]
    rotation_tensor = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensor = zeros(2, 2)
    tensor[1, 1] = 1
    dof = 2
    back = true
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, back)
    @test tensorTest == tensor
    angles = [90.0]
    rotation_tensor = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, back)
    @test isapprox(tensorTest[1, 1] + 1, 1) # plus one, because of how approx works
    @test isapprox(tensorTest[1, 2] + 1, 1)
    @test isapprox(tensorTest[2, 1] + 1, 1)
    # @test isapprox(tensorTest[2, 2], 1)
    back = false
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, back)
    @test tensorTest == tensor

    angles = [0, 0, 0]
    rotation_tensor = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensor = zeros(3, 3)
    tensor[1, 1] = 1
    dof = 3

    back = true
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, back)
    @test tensorTest == tensor
    angles = [0, 0, 90.0]
    rotation_tensor = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, back)
    @test isapprox(tensorTest[1, 1] + 1, 1) # plus one, because of how approx works
    @test isapprox(tensorTest[1, 2] + 1, 1)
    @test isapprox(tensorTest[1, 3] + 1, 1)
    @test isapprox(tensorTest[2, 1] + 1, 1)
    # @test isapprox(tensorTest[2, 2], 1)
    @test isapprox(tensorTest[2, 3] + 1, 1)
    @test isapprox(tensorTest[3, 1] + 1, 1)
    @test isapprox(tensorTest[3, 2] + 1, 1)
    @test isapprox(tensorTest[3, 3] + 1, 1)

    back = false
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, back)
    @test tensorTest == tensor

    angles = [10, 20, 90.0]
    rotation_tensor = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, true)
    tensorTest = Correspondence.rotate_second_order_tensor(rotation_tensor, tensor, dof, false)
    @test tensorTest == tensor

end