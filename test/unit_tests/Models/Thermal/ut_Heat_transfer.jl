# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

#include("../../../../src/PeriLab.jl")
#import .PeriLab

@test parentmodule(PeriLab.ParameterSpec.lookup_model(:thermal, "Heat Transfer")) ===
      PeriLab.Solver_Manager.Model_Factory.Thermal.Heat_Transfer

HT = PeriLab.Solver_Manager.Model_Factory.Thermal.Heat_Transfer

@testset "ut_specific_volume_grid_2d" begin
    # 3x3 grid, spacing 1, all nodes are neighbors of each other
    coor = [Float64(i % 3) for i in 0:8]
    coor = hcat(coor, [Float64(i ÷ 3) for i in 0:8])
    nodes = collect(1:9)
    nlist = [[j for j in 1:9 if j != i] for i in 1:9]
    bond_geometry = [[coor[j, :] - coor[i, :] for j in nlist[i]] for i in 1:9]
    extent = ones(9, 2)
    bond_face = [zeros(Int64, 8) for _ in 1:9]
    bond_face_area = [zeros(8) for _ in 1:9]
    HT.init_bond_faces!(bond_face, bond_face_area, nodes, nlist, bond_geometry, extent,
                        nothing, Val(2))
    # diagonal and second neighbors cover nothing
    @test count(!=(0), bond_face[5]) == 4
    @test count(!=(0), bond_face[1]) == 2
    bed_face = zeros(Int64, 9)
    HT.init_bed_faces!(bed_face, nodes, coor, extent, nothing, -0.5, Val(2))
    @test bed_face == [4, 4, 4, 0, 0, 0, 0, 0, 0]

    sv = zeros(9)
    bed_sv = zeros(9)
    surface = fill(true, 9)
    active = fill(true, 9)
    damage = [ones(8) for _ in 1:9]
    HT.compute_specific_volume!(sv, bed_sv, surface, nodes, nlist, active, damage,
                                bond_face, bond_face_area, extent, bed_face, true, Val(2))
    @test sv == [1, 0, 1, 1, 0, 1, 2, 1, 2]
    @test bed_sv == [1, 1, 1, 0, 0, 0, 0, 0, 0]
    @test surface == [true, true, true, true, false, true, true, true, true]

    # an inactive node and a broken bond open the faces again
    active[2] = false
    damage[5][findfirst(==(6), nlist[5])] = 0
    surface .= true
    HT.compute_specific_volume!(sv, bed_sv, surface, nodes, nlist, active, damage,
                                bond_face, bond_face_area, extent, bed_face, false, Val(2))
    @test sv[5] == 2
    @test sv[1] == 2
    @test bed_sv[1] == 1

    heat_flow = zeros(9)
    temperature = fill(10.0, 9)
    # large heat capacity: the relaxation factor is one
    HT.apply_surface_heat_flow!(heat_flow, temperature, sv, bed_sv, surface, nodes, 2.0,
                                0.0, 3.0, 4.0, ones(9), fill(1e12, 9), 1.0)
    @test heat_flow[1] ≈ 2.0 * 10 * 2 + 3.0 * 6 * 1
    @test heat_flow[5] ≈ 2.0 * 10 * 2
end

@testset "ut_surface_heat_flow_exact_relaxation" begin
    # a = 2 * 1 + 3 * 1 = 5 towards Tw = (2 * 0 + 3 * 4) / 5; x = 5 * 10 / 1 is far beyond
    # the explicit limit 2, the explicit update must still give the exact solution
    heat_flow = [0.0]
    T = 10.0
    HT.apply_surface_heat_flow!(heat_flow, [T], [1.0], [1.0], [true], [1], 2.0, 0.0, 3.0,
                                4.0, [1.0], [1.0], 10.0)
    Tw = 12 / 5
    @test T - heat_flow[1] * 10.0 / 1.0 ≈ Tw + (T - Tw) * exp(-50)
    @test HT.relaxation_factor(0.0) == 1
    @test HT.relaxation_factor(1e-12) ≈ 1
end

@testset "ut_specific_volume_crossing_strands_3d" begin
    # G-code strands (sampling 2, width 0.42, height 0.2); the upper one is rotated by 90°
    extent = [2.0 0.42 0.2; 2.0 0.42 0.2]
    rotation = zeros(2, 3, 3)
    rotation[1, :, :] = [1 0 0; 0 1 0; 0 0 1]
    rotation[2, :, :] = [0 -1 0; 1 0 0; 0 0 1]
    nlist = [[2], [1]]
    bond_geometry = [[[0.0, 0.0, 0.2]], [[0.0, 0.0, -0.2]]]
    bond_face = [[0], [0]]
    bond_face_area = [[0.0], [0.0]]
    HT.init_bond_faces!(bond_face, bond_face_area, [1, 2], nlist, bond_geometry, extent,
                        rotation, Val(3))
    @test bond_face == [[5], [6]]
    @test bond_face_area[1][1] ≈ 0.42^2
    @test bond_face_area[2][1] ≈ 0.42^2

    bed_face = [0, 0]
    HT.init_bed_faces!(bed_face, [1, 2], [0.0 0.0 0.1; 0.0 0.0 0.3], extent, rotation,
                       0.0, Val(3))
    @test bed_face == [6, 0]

    sv = zeros(2)
    bed_sv = zeros(2)
    HT.compute_specific_volume!(sv, bed_sv, [true, true], [1, 2], nlist, [true, true],
                                [[1.0], [1.0]], bond_face, bond_face_area, extent,
                                bed_face, true, Val(3))
    box = 2 * 0.42 * 0.2
    sides = 2 * (0.42 * 0.2 + 2 * 0.2)
    top = 2 * 0.42
    @test sv[1] ≈ (sides + top - 0.42^2) / box
    @test bed_sv[1] ≈ top / box
    @test sv[2] ≈ (sides + 2top - 0.42^2) / box
    @test bed_sv[2] == 0
end

@testset "ut_surface_update_nodes" begin
    # chain 1-2-3-4, node 5 is a ghost neighbor of node 4
    nlist = [[2], [1, 3], [2, 4], [3, 5], [4]]
    buffer = zeros(Int64, 5)
    update = fill(false, 5)
    surface_active = fill(false, 5)
    border = [false, false, false, true, false]
    active = [true, true, false, false, false]
    @test HT.surface_update_nodes!(buffer, update, surface_active, border, active, [1, 2],
                                   nlist, 4) == [1, 2]
    # nothing changed: nothing to update
    @test isempty(HT.surface_update_nodes!(buffer, update, surface_active, border, active,
                                           [1, 2], nlist, 4))
    # node 3 is activated: its neighbors 2 and 4 change, 4 only once it is active
    active[3] = true
    @test HT.surface_update_nodes!(buffer, update, surface_active, border, active,
                                   [1, 2, 3], nlist, 4) == [2, 3]
    active[4] = true
    @test HT.surface_update_nodes!(buffer, update, surface_active, border, active,
                                   [1, 2, 3, 4], nlist, 4) == [3, 4]
    # an activated ghost node updates the border nodes
    active[5] = true
    @test HT.surface_update_nodes!(buffer, update, surface_active, border, active,
                                   [1, 2, 3, 4], nlist, 4) == [4]
end

@testset "ut_compute_thermal_model" begin
    # three nodes in a row on the bed at y = -0.5
    DM = PeriLab.Data_Manager
    DM.initialize_data()
    DM.set_num_controller(3)
    DM.set_dof(2)
    nn = DM.create_constant_node_scalar_field("Number of Neighbors", Int64)
    nn .= [1, 2, 1]
    nlist = DM.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [2]
    coor = DM.create_constant_node_vector_field("Coordinates", Float64, 2)
    coor[:, 1] = [0.0, 1.0, 2.0]
    bond_geometry = DM.create_constant_bond_vector_state("Bond Geometry", Float64, 2)
    for iID in 1:3, (jID, neighborID) in enumerate(nlist[iID])
        bond_geometry[iID][jID] = coor[neighborID, :] - coor[iID, :]
    end
    DM.create_constant_node_scalar_field("Volume", Float64; default_value = 1)
    # large heat capacity: the relaxation factor of the surface heat flow is one
    DM.create_constant_node_scalar_field("Density", Float64; default_value = 1)
    DM.create_constant_node_scalar_field("Specific Heat Capacity", Float64;
                                         default_value = 1e12)
    DM.create_constant_node_scalar_field("Active", Bool; default_value = true)
    DM.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    DM.create_node_scalar_field("Heat Flow", Float64)
    DM.create_node_scalar_field("Temperature", Float64)
    DM.get_field("Temperature", "NP1") .= [0, 10, 0]
    DM.create_constant_node_scalar_field("Specific Volume", Float64)
    DM.create_constant_node_scalar_field("Surface_Nodes", Bool; default_value = true)
    th = typed_model(:thermal,
                     Dict("Thermal Model" => "Heat Transfer",
                          "Heat Transfer Coefficient" => 2,
                          "Environmental Temperature" => 0,
                          "Bed Temperature" => "4",
                          "Bed Position" => -0.5);
                     name_key = "Thermal Model")
    HT.init_model(collect(1:3), th.model, th, 1)
    HT.compute_model(collect(1:3), th.model, th, 1, 0.0, 1.0)
    # end nodes: outer end and top exposed, middle node: top only
    @test DM.get_field("Specific Volume") == [2, 1, 2]
    @test DM.get_field("Bed Specific Volume") == [1, 1, 1]
    @test DM.get_field("Heat Flow", "NP1") ≈ [2 * (0 - 4), 2 * 10 + 2 * (10 - 4),
        2 * (0 - 4)]

    # conduction to the bed across half the node height (1): coefficient 2 * 1.5 / 1
    th = typed_model(:thermal,
                     Dict("Thermal Model" => "Heat Transfer",
                          "Heat Transfer Coefficient" => 2,
                          "Environmental Temperature" => 0,
                          "Bed Temperature" => 4,
                          "Bed Thermal Conductivity" => 1.5,
                          "Bed Position" => -0.5);
                     name_key = "Thermal Model")
    DM.get_field("Heat Flow", "NP1") .= 0
    HT.compute_model(collect(1:3), th.model, th, 1, 0.0, 1.0)
    @test DM.get_field("Heat Flow", "NP1") ≈ [3 * (0 - 4), 2 * 10 + 3 * (10 - 4),
        3 * (0 - 4)]
end
