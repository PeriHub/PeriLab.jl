# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
#using .PeriLab

@testset "material_name" begin
    @test PeriLab.Solver_Manager.Model_Factory.Material.Unified_Bondbased_Elastic.material_name() ==
          "Unified Bond-based Elastic"
    @test !(PeriLab.Solver_Manager.Model_Factory.Material.Unified_Bondbased_Elastic.fe_support())
end
@testset "compute_bond_based_strain" begin

    # Initialize input arrays and parameters
    bb_strain = zeros(3, 3)
    deformed_bond = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    undeformed_bond = [1.0 1.5 2.0; 2.5 3.0 3.5; 4.0 4.5 5.0]
    bond_damage = [0.1, 0.2, 0.3]
    nlist = [1, 2, 3]
    volume = [1.0, 1.0, 1.0]

    # # Call the function
    # PeriLab.Solver_Manager.Model_Factory.Material.Unified_Bondbased_Elastic.compute_bond_based_strain(bb_strain,
    #                                                                                                   deformed_bond,
    #                                                                                                   undeformed_bond,
    #                                                                                                   bond_damage,
    #                                                                                                   nlist,
    #                                                                                                   volume)

    # # Assert expected results
    # @test bb_strain == [0.575 0.575 0.575
    #                     0.575 0.575 0.575
    #                     0.575 0.575 0.575]
end
