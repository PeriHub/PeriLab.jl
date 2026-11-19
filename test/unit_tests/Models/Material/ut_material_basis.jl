# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
#include("../../../../src/PeriLab.jl")
#using .PeriLab

@testset "ut_init_local_damping_due_to_damage" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.create_constant_node_scalar_field("Horizon", Float64) .= 1.0
    local_damping = (representative_youngs_modulus = 70.0e9, damping_coefficient = 0.5)
    PeriLab.Solver_Manager.Material_Basis.init_local_damping_due_to_damage(collect(1:2),
                                                                           "plane strain",
                                                                           local_damping)
    @test PeriLab.Data_Manager.has_key("Bond Based Constant")
end

@testset "ut_local_damping_due_to_damage" begin
    # two nodes on the x axis, node 1 damaged, bond elongating
    DM = PeriLab.Data_Manager
    DM.initialize_data()
    DM.set_num_controller(2)
    DM.set_dof(2)
    DM.create_constant_node_scalar_field("Number of Neighbors", Int64) .= 1
    nlist = DM.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] .= 2
    nlist[2] .= 1
    DM.create_constant_node_scalar_field("Horizon", Float64) .= 3.0
    DM.create_constant_node_scalar_field("Density", Float64) .= 1.0
    DM.create_constant_node_scalar_field("Volume", Float64) .= [1.0, 2.0]
    _, damage = DM.create_node_scalar_field("Damage", Float64)
    damage .= [0.5, 0.0]
    lN, lNP1 = DM.create_bond_scalar_state("Deformed Bond Length", Float64)
    for i in 1:2
        lN[i] .= 1.0
        lNP1[i] .= 1.1
    end
    _, geo = DM.create_bond_vector_state("Deformed Bond Geometry", Float64, 2)
    geo[1][1] = [1.1, 0.0]
    geo[2][1] = [-1.1, 0.0]
    _, force = DM.create_node_vector_field("Force Densities", Float64, 2)
    params = (representative_youngs_modulus = 4.0, damping_coefficient = 0.5)
    MB = PeriLab.Solver_Manager.Material_Basis
    MB.init_local_damping_due_to_damage([1, 2], "plane stress", params)
    MB.local_damping_due_to_damage([1, 2], params, 0.1)

    c = MB.bond_based_constant("plane stress", 3.0)
    # each bond visit pulls the nodes together with c*E*k*avg_damage*rate/v0
    t = 0.5 * 4.0 * c * 0.25 * (0.1 / 0.1) / 2.0
    @test force[1, 1] ≈ 2t * 2.0 # pulls node 1 towards node 2 (opposes elongation)
    @test force[2, 1] ≈ -2t * 1.0
    @test force[1, 1] * 1.0 + force[2, 1] * 2.0 ≈ 0 atol=1e-14 # momentum
    @test all(force[:, 2] .== 0)

    # zero coefficient changes nothing
    force .= 0
    MB.local_damping_due_to_damage([1, 2], (; params..., damping_coefficient = 0.0), 0.1)
    @test all(force .== 0)

    # energy: broken bond, nodes flying apart, damping is the only force
    damage .= 1.0
    DM.get_field("Volume") .= 1.0
    x = [0.0, 1.0]
    v = [-1.0, 1.0]
    dt = 0.1
    energy = [sum(v .^ 2) / 2]
    for _ in 1:50
        lN[1][1] = lN[2][1] = x[2] - x[1]
        x .+= v .* dt
        lNP1[1][1] = lNP1[2][1] = x[2] - x[1]
        geo[1][1] = [x[2] - x[1], 0.0]
        geo[2][1] = [x[1] - x[2], 0.0]
        force .= 0
        MB.local_damping_due_to_damage([1, 2], params, dt)
        v .+= force[:, 1] .* dt # density 1
        push!(energy, sum(v .^ 2) / 2)
    end
    @test all(diff(energy) .< 0)
    @test energy[end] < 0.5 * energy[1]
    @test sum(v) ≈ 0 atol=1e-14

    # no damping without damage
    force .= 0
    damage .= 0
    MB.local_damping_due_to_damage([1, 2], params, 0.1)
    @test all(force .== 0)
end

@testset "ut_apply_pointwise_E" begin
    nodes = 2:3

    bond_force = [[ones(2), ones(2)], [ones(2), ones(2)], [ones(2), ones(2)]]
    E = 3.3
    bond_force[3][2][2] = 3
    bond_force = bond_field(bond_force)
    PeriLab.Solver_Manager.Material_Basis.apply_pointwise_E(nodes, E, bond_force)
    @test bond_force[1][1][1] == 1
    @test bond_force[1][1][2] == 1
    @test bond_force[1][2][1] == 1
    @test bond_force[1][2][2] == 1
    @test bond_force[2][2][1] == E
    @test bond_force[2][2][2] == E
    @test bond_force[2][2][1] == E
    @test bond_force[2][2][2] == E
    @test bond_force[3][1][1] == E
    @test bond_force[3][1][2] == E
    @test bond_force[3][2][1] == E
    @test bond_force[3][2][2] == 3 * E

    bond_force = bond_field([[ones(2), ones(2)], [ones(2), ones(2)],
                                [ones(2), ones(2), ones(2)]])
    E = zeros(4)
    PeriLab.Solver_Manager.Material_Basis.apply_pointwise_E(nodes, E, bond_force)
    @test bond_force[1][1][1] == 1
    @test bond_force[1][1][2] == 1
    @test bond_force[1][2][1] == 1
    @test bond_force[1][2][2] == 1
    @test bond_force[2][1][1] == 0
    @test bond_force[2][1][2] == 0
    @test bond_force[2][2][1] == 0
    @test bond_force[2][2][2] == 0
    @test bond_force[3][1][1] == 0
    @test bond_force[3][1][2] == 0
    @test bond_force[3][2][1] == 0
    @test bond_force[3][2][2] == 0
    @test bond_force[3][3][1] == 0
    @test bond_force[3][3][2] == 0
end

@testset "ut_distribute_forces" begin
    nodes = [1, 2]
    dof = 2
    nlist = [[2, 3], [1, 3]]
    nBonds = fill(dof, 2)
    bond_force = [[fill(1.0, dof) for j in 1:n] for n in nBonds]
    volume = [1.0, 1.0, 1.0]
    bond_damage = [fill(0.5, n) for n in nBonds]
    force_densities = zeros(3, 2)

    expected_force_densities = copy(force_densities)
    for iID in nodes
        expected_force_densities[iID,
                                 :] .+= transpose(sum(bond_damage[iID] .*
                                                      mapreduce(permutedims, vcat,
                                                                bond_force[iID]) .*
                                                      volume[nlist[iID]],
                                                      dims = 1))
        expected_force_densities[nlist[iID],
                                 :] .-= bond_damage[iID] .*
                                        mapreduce(permutedims, vcat,
                                                  bond_force[iID]) .*
                                        volume[iID]
    end

    PeriLab.Solver_Manager.Material_Basis.distribute_forces!(force_densities, nodes,
                                                             bond_field(nlist),
                                                             bond_field(bond_force), volume,
                                                             bond_field(bond_damage))
    @test force_densities ≈ expected_force_densities
end
@testset "ut_compute_Piola_Kirchhoff_stress" begin
    stress = [1.0 0.0; 0.0 1.0]
    deformation_gradient = [2.0 0.0; 0.0 2.0]
    expected_result = [2.0 0.0; 0.0 2.0]
    result = zeros(2, 2)
    PeriLab.Solver_Manager.Material_Basis.compute_Piola_Kirchhoff_stress!(result, stress,
                                                                          deformation_gradient)
    @test isapprox(result, expected_result)
end

@testset "ut_local_damping_first_step" begin
    # damaged from t = 0 (e.g. bond filter with contact), nothing moves: no damping force
    DM = PeriLab.Data_Manager
    DM.initialize_data()
    DM.set_num_controller(3)
    DM.set_dof(2)
    DM.create_constant_node_scalar_field("Number of Neighbors", Int64) .= 2
    nlist = DM.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] .= [2, 3]
    nlist[2] .= [1, 3]
    nlist[3] .= [1, 2]
    x = [0.0 0.0; 1.0 0.0; 0.0 1.0]
    geometry = DM.create_constant_bond_vector_state("Bond Geometry", Float64, 2)
    len = DM.create_constant_bond_scalar_state("Bond Length", Float64)
    for i in 1:3, (j, n) in enumerate(nlist[i])
        geometry[i][j] = x[n, :] - x[i, :]
        len[i][j] = hypot((x[n, :] - x[i, :])...)
    end
    DM.create_deformed_bond_states(2)
    for time in ("N", "NP1")
        @test DM.get_field("Deformed Bond Length", time) == len
        @test DM.get_field("Deformed Bond Geometry", time) == geometry
    end
    DM.create_constant_node_scalar_field("Horizon", Float64) .= 1.5
    DM.create_constant_node_scalar_field("Density", Float64) .= 2.0
    DM.create_constant_node_scalar_field("Volume", Float64) .= 0.5
    _, damage = DM.create_node_scalar_field("Damage", Float64)
    damage .= [0.4, 0.1, 0.0]
    _, force = DM.create_node_vector_field("Force Densities", Float64, 2)
    params = (representative_youngs_modulus = 8.0, damping_coefficient = 0.3)
    MB = PeriLab.Solver_Manager.Material_Basis
    MB.init_local_damping_due_to_damage(1:3, "plane stress", params)
    MB.local_damping_due_to_damage(1:3, params, 0.1)
    @test all(force .== 0)

    # a bond without deformed state is skipped, not 0/0
    DM.get_field("Deformed Bond Length", "N")[1][1] = 0.0
    MB.local_damping_due_to_damage(1:3, params, 0.1)
    @test all(isfinite, force)

    # stability number c E k ΣV dt / (v0 ρ)
    k = MB.bond_based_constant("plane stress", 1.5)
    @test MB.local_damping_stability(1:3, params, 0.1) ≈
          0.3 * 8.0 * k * 1.0 * 0.1 / (2.0 * 2.0)
end
