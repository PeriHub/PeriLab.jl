# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test
include(
    "../../../../../../src/Models/Material/Material_Models/Ordinary/PD_Solid_Plastic.jl",
)
#include("../../../../../../src/PeriLab.jl")
#using .PeriLab

@testset "get_name&fe_support" begin
    @test PD_Solid_Plastic.material_name() == "PD Solid Plastic"
    @test !(PD_Solid_Plastic.fe_support())
end

@testset "ut_init_material_model" begin

    nodes = 2
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nodes)
    dof = 3
    test_data_manager.set_dof(dof)
    horizon = test_data_manager.create_constant_node_field("Horizon", Float64, 1)

    horizon[1] = 3
    horizon[2] = 2

    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)

    nn .= 1
    @test isnothing(
        PD_Solid_Plastic.init_material_model(
            test_data_manager,
            Vector{Int64}(1:nodes),
            Dict(),
        ),
    )

    test_data_manager = PD_Solid_Plastic.init_material_model(
        test_data_manager,
        Vector{Int64}(1:nodes),
        Dict("Yield Stress" => 5.3),
    )
    yield = test_data_manager.get_field("Yield Value")

    @test isapprox(yield[1], 25 * 5.3 * 5.3 / (8 * pi * 3^5))
    @test isapprox(yield[2], 25 * 5.3 * 5.3 / (8 * pi * 2^5))

    test_data_manager = PD_Solid_Plastic.init_material_model(
        test_data_manager,
        Vector{Int64}(1:nodes),
        Dict("Yield Stress" => 2.2, "Symmetry" => "plane stress"),
    )
    yield = test_data_manager.get_field("Yield Value")

    @test isapprox(yield[1], 225 * 2.2 * 2.2 / (24 * pi * 3^4))
    @test isapprox(yield[2], 225 * 2.2 * 2.2 / (24 * pi * 2^4))

end
