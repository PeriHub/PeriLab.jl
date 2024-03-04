# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../src/Physics/Additive/Additive_Factory.jl")
include("../../../../src/Support/data_manager.jl")


using Test
using .Additive
@testset "init_additive" begin
    test_Data_manager = Data_manager
    test_Data_manager.properties[23] = Dict("Additive Model" => Dict("Additive Model" => "does not exist"))
    println()
    @test isnothing(Additive.init_additive_model(test_Data_manager, Vector{Int64}([1, 2, 3]), 23))
end