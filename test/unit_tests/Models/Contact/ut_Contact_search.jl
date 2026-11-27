# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
using Random: Xoshiro

const UT_CS = PeriLab.Solver_Manager.Model_Factory.Contact.Contact_Search

# sorted (master, slave) pairs within `radius`, brute force
function ut_cs_brute_pairs(positions, masters, slaves, radius)
    found = Tuple{Int64,Int64}[]
    for m in masters, s in slaves
        sum(abs2, positions[m, :] .- positions[s, :]) <= radius^2 && push!(found, (m, s))
    end
    return sort!(found)
end

# sorted (master, slave) pairs of a CSR list
function ut_cs_pairs(masters, offsets, slaves)
    found = Tuple{Int64,Int64}[]
    for (k, m) in enumerate(masters), c in offsets[k]:(offsets[k + 1] - 1)
        push!(found, (m, slaves[c]))
    end
    return sort!(found)
end
ut_cs_pairs(p) = ut_cs_pairs(p.masters, p.offsets, p.slaves)

function ut_cs_setup(positions)
    PeriLab.Data_Manager.initialize_data()
    n = size(positions, 1)
    PeriLab.Data_Manager.set_dof(size(positions, 2))
    PeriLab.Data_Manager.set_num_controller(n)
    PeriLab.Data_Manager.create_constant_node_scalar_field("Contact Nodes", Int64)
    PeriLab.Data_Manager.set_exchange_id_to_local_id(Dict(i => i for i in 1:n))
    PeriLab.Data_Manager.set_all_positions(positions)
end

@testset "ut_contact_search_brute_force" begin
    rng = Xoshiro(42)
    for dof in (2, 3)
        positions = rand(rng, 60, dof)
        masters = collect(1:30)
        slaves = collect(31:60)
        ut_cs_setup(positions)
        s = UT_CS.ContactSearch(dof, masters, slaves, 0.2, 0.3, 1)
        UT_CS.compute_contact_pairs!(s)
        p = s.pairs
        expected = ut_cs_brute_pairs(positions, masters, slaves, 0.2)
        @test !isempty(expected)
        @test ut_cs_pairs(p) == expected
        # candidates: all pairs within the search radius
        @test ut_cs_pairs(s.candidate_masters, s.candidate_offsets, s.candidate_slaves) ==
              ut_cs_brute_pairs(positions, masters, slaves, 0.3)
        @test s.candidate_nodes == sort!(unique(vcat(s.candidate_masters,
                                s.candidate_slaves)))
        # distances and unit normals along the master-slave direction
        for (k, m) in enumerate(p.masters), c in p.offsets[k]:(p.offsets[k + 1] - 1)
            diff = positions[p.slaves[c], :] .- positions[m, :]
            normal = p.normals[((c - 1) * dof + 1):(c * dof)]
            @test p.distances[c] ≈ sqrt(sum(abs2, diff))
            @test abs(sum(normal .* diff)) ≈ p.distances[c]
            @test sum(abs2, normal) ≈ 1
        end
        contact_nodes = PeriLab.Data_Manager.get_field("Contact Nodes")
        @test all(contact_nodes[p.masters] .== 5)
        @test all(contact_nodes[p.slaves] .== 8)
    end
end

@testset "ut_contact_search_frequency" begin
    # master at the origin, slaves on the x axis
    positions = [0.0 0.0; 0.15 0.0; 0.25 0.0; 0.5 0.0]
    ut_cs_setup(positions)
    s = UT_CS.ContactSearch(2, [1], [2, 3, 4], 0.2, 0.3, 3)
    @test UT_CS.global_search_due(s)
    UT_CS.compute_contact_pairs!(s)
    @test ut_cs_pairs(s.pairs) == [(1, 2)]
    @test s.candidate_slaves == [2, 3]
    # node 3 (candidate) comes into contact: found without a global search;
    # node 4 (no candidate) is missed until the next global search
    positions[3, 1] = 0.1
    positions[4, 1] = 0.05
    @test !UT_CS.global_search_due(s)
    UT_CS.compute_contact_pairs!(s)
    @test ut_cs_pairs(s.pairs) == [(1, 2), (1, 3)]
    UT_CS.compute_contact_pairs!(s)
    @test ut_cs_pairs(s.pairs) == [(1, 2), (1, 3)]
    @test UT_CS.global_search_due(s)
    UT_CS.compute_contact_pairs!(s)
    @test ut_cs_pairs(s.pairs) == [(1, 2), (1, 3), (1, 4)]
end

@testset "ut_contact_search_empty_and_radius" begin
    ut_cs_setup([0.0 0.0; 1.0 0.0])
    s = UT_CS.ContactSearch(2, Int64[], [2], 0.2, 0.3, 1)
    UT_CS.compute_contact_pairs!(s)
    @test isempty(s.pairs.masters)
    @test s.pairs.offsets == [1]
    # a search radius below the contact radius is raised to it
    @test UT_CS.ContactSearch(2, [1], [2], 0.2, 0.1, 1).search_radius == 0.2
    # runtests.jl disables warnings for the unit tests
    Logging.disable_logging(Logging.Info)
    try
        @test_logs (:warn, r"Search Radius") UT_CS.ContactSearch(2, [1], [2], 0.2, 0.2, 2)
        @test_logs UT_CS.ContactSearch(2, [1], [2], 0.2, 0.3, 2)
        @test_logs UT_CS.ContactSearch(2, [1], [2], 0.2, 0.2, 1)
    finally
        Logging.disable_logging(Logging.Warn)
    end
end
