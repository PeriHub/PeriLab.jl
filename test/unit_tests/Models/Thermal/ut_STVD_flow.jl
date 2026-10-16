# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const STVD = PeriLab.Solver_Manager.Model_Factory.Thermal.STVD_Flow

@test STVD.thermal_model_name() == "STVD Flow"

@testset "ut_compute_depth" begin
    # chain of 7 nodes, spacing 1, nodes 3:5 are diffusant (block 2)
    coordinates = zeros(7, 2)
    coordinates[:, 1] .= 1:7
    block_ids = [1, 1, 2, 2, 2, 1, 1]
    bond_length = [[1.0] for _ in 1:7]
    depth = zeros(7)
    STVD.compute_depth!(depth, coordinates, block_ids, [2], bond_length)
    @test depth ≈ [-1, -1, 0.5, 1.5, 0.5, -1, -1]
    @test STVD.parse_block_ids("2 3") == [2, 3]
    @test STVD.parse_block_ids(2) == [2]
end

@testset "ut_compute_stvd_flow" begin
    nlist = [[2], [1, 3], [2, 4], [3, 5], [4]]
    bond_length = [ones(length(n)) for n in nlist]
    bond_damage = [ones(length(n)) for n in nlist]
    depth = [-1.0, 0.5, 1.5, 0.5, -1.0]
    concentration = [0.0, 1.0, 1.0, 1.0, 0.0]
    volume = ones(5)
    horizon = ones(5)
    spacing = ones(5)

    # front not yet in the core -> only surface bonds of the core are active
    flow = zeros(5)
    STVD.compute_stvd_flow!(1:5, 2, nlist, 1.0, 0.6, false, depth, spacing, bond_damage,
                            bond_length, horizon, concentration, volume, flow)
    @test flow[3] == 0
    @test flow[2] ≈ 6 / pi
    @test flow[1] ≈ -6 / pi
    @test isapprox(sum(flow), 0; atol = 1e-14)

    # no reaction -> core is frozen, substrate has no gradient
    flow .= 0
    STVD.compute_stvd_flow!(1:5, 2, nlist, 1.0, 0.0, false, depth, spacing, bond_damage,
                            bond_length, horizon, concentration, volume, flow)
    @test all(flow .== 0)
end

@testset "ut_compute_front_depth" begin
    @test STVD.compute_front_depth(Dict("Reaction Rate" => 0.1), 2.5) ≈ 0.25
    @test STVD.compute_front_depth(Dict("Reaction Rate" => 0.1, "Reaction Law" => "Linear"),
                                   2.5) ≈ 0.25
    parameter = Dict("Reaction Rate" => 0.25, "Reaction Law" => "Exponential",
                     "Core Size" => 0.5)
    @test STVD.compute_front_depth(parameter, 0.0) == 0
    @test STVD.compute_front_depth(parameter, 4.0) ≈ 0.5 * (1 - exp(-1))
end

@testset "ut_init_model_errors" begin
    PeriLab.Data_Manager.initialize_data()
    nodes = Vector{Int64}(1:3)
    @test_logs (:error, "Thermal Conductivity not defined.") @test_throws PeriLab.PeriLabError STVD.init_model(nodes,
                                                                                                               Dict{String,
                                                                                                                    Any}())
    @test_logs (:error, "Reaction Rate not defined.") @test_throws PeriLab.PeriLabError STVD.init_model(nodes,
                                                                                                        Dict{String,
                                                                                                             Any}("Thermal Conductivity" =>
                                                                                                                      1.0))
    parameter = Dict{String,Any}("Thermal Conductivity" => 1.0, "Reaction Rate" => 0.1,
                                 "Diffusant Block ID" => 2, "Reaction Law" => "Cubic")
    @test_logs (:error,
                "Reaction Law Cubic not supported. Use ''Linear'' or ''Exponential''.") @test_throws PeriLab.PeriLabError STVD.init_model(nodes,
                                                                                                                                          parameter)
    parameter["Reaction Law"] = "Exponential"
    @test_logs (:error, "Core Size has to be defined for the exponential reaction law.") @test_throws PeriLab.PeriLabError STVD.init_model(nodes,
                                                                                                                                           parameter)
    PeriLab.Data_Manager.data["mpi_active"] = true
    @test_logs (:error,
                "STVD Flow does not support MPI yet, the depth below the diffusant surface needs all substrate nodes.") @test_throws PeriLab.PeriLabError STVD.init_model(nodes,
                                                                                                                                                                          parameter)
    PeriLab.Data_Manager.data["mpi_active"] = false
end

@testset "ut_init_model_depth_from_mesh" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    depth = PeriLab.Data_Manager.create_constant_node_scalar_field("STVD_Depth", Float64)
    depth .= [-1.0, 0.25, 0.5]
    STVD.init_model(Vector{Int64}(1:3),
                    Dict{String,Any}("Thermal Conductivity" => 1.0, "Reaction Rate" => 0.1,
                                     "Diffusant Block ID" => 2))
    @test PeriLab.Data_Manager.get_field("STVD_Depth") == [-1.0, 0.25, 0.5]
end

@testset "ut_continuous_release" begin
    # node-wise release
    @test STVD.inert_fraction(0.5, 1.0, 0.4, false) == 1
    @test STVD.inert_fraction(0.5, 1.0, 0.5, false) == 0
    # continuous release: depth 0.5, spacing 1 -> node covers depths 0 ... 1
    @test STVD.inert_fraction(0.5, 1.0, 0.0, true) == 1
    @test STVD.inert_fraction(0.5, 1.0, 0.25, true) ≈ 0.75
    @test STVD.inert_fraction(0.5, 1.0, 1.0, true) == 0
    @test STVD.inert_fraction(-1.0, 1.0, 0.0, true) == 0

    # chain substrate - core - substrate, front inside the outer core nodes
    nlist = [[2], [1, 3], [2, 4], [3, 5], [4]]
    bond_length = [ones(length(n)) for n in nlist]
    bond_damage = [ones(length(n)) for n in nlist]
    depth = [-1.0, 0.5, 1.5, 0.5, -1.0]
    spacing = ones(5)
    concentration = [0.0, 1.0, 1.0, 1.0, 0.0]
    volume = ones(5)
    horizon = ones(5)
    flow = zeros(5)
    STVD.compute_stvd_flow!(1:5, 2, nlist, 1.0, 0.25, true, depth, spacing, bond_damage,
                            bond_length, horizon, concentration, volume, flow)
    # outer core nodes are 25 % released: mobile 0.25 flows into the substrate
    @test flow[1] ≈ -6 / pi * 0.25
    @test flow[2] ≈ 6 / pi * 0.25
    @test flow[3] == 0
    @test isapprox(sum(flow), 0; atol = 1e-14)
end
