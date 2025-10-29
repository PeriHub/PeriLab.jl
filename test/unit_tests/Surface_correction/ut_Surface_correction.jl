# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

#include("../../../src/PeriLab.jl")
#import .PeriLab

@testset "ut_init_surface_correction" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_block_id_list([1])
    test_data_manager.init_properties()
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(4)
    block_iD = test_data_manager.create_constant_node_scalar_field("Block_Id", Int64)
    block_iD .= 1
    mod_struct = PeriLab.Solver_Manager.Model_Factory

    mod_struct.init_surface_correction(Dict(),
                                       "local_synch",
                                       "synchronise_field")
    @test isnothing(mod_struct.init_surface_correction(Dict("Surface Correction" => Dict("a" => 0)),
                                                       "local_synch",
                                                       "synchronise_field"))
end
