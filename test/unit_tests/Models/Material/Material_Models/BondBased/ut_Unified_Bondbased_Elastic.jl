# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include(
    "../../../../../../src/Models/Material/Material_Models/BondBased/Unified_Bondbased_Elastic.jl",
)

using .Unified_Bondbased_Elastic
using Test
using TimerOutputs
#include("../../../../../../src/PeriLab.jl")
#using .PeriLab

# const to = TimerOutput()

@testset "material_name" begin
    @test Unified_Bondbased_Elastic.material_name() == "Unified Bond-based Elastic"
    @test !(Unified_Bondbased_Elastic.fe_support())
end
@testset "compute_bond_based_strain" begin

    # Initialize input arrays and parameters
    bb_strain = zeros(3, 3)
    deformed_bond = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    undeformed_bond = [1.0 1.5 2.0; 2.5 3.0 3.5; 4.0 4.5 5.0]
    bond_damage = [0.1, 0.2, 0.3]
    nlist = [1, 2, 3]
    volume = [1.0, 1.0, 1.0]

    # Call the function
    Unified_Bondbased_Elastic.compute_bond_based_strain(
        bb_strain,
        deformed_bond,
        undeformed_bond,
        bond_damage,
        nlist,
        volume,
    )

    # Assert expected results
    @test bb_strain == [
        0.575 0.4999999999999999 0.4428571428571428
        0.7874999999999999 0.6666666666666667 0.5821428571428571
        0.9999999999999998 0.8333333333333333 0.7214285714285713
    ]


end
