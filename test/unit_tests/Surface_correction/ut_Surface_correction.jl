# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

#include("../../../src/PeriLab.jl")
#import .PeriLab

@testset "ut_init_surface_correction" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_block_id_list([1])
    PeriLab.Data_Manager.init_properties()
    PeriLab.Data_Manager.set_dof(3)
    PeriLab.Data_Manager.set_num_controller(4)
    block_iD = PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64)
    block_iD .= 1
    mod_struct = PeriLab.Solver_Manager.Model_Factory

    mod_struct.init_surface_correction(Dict(),
                                       "local_synch",
                                       "synchronise_field")
    @test_logs (:error, "Surface Correction needs a Type definition") mod_struct.init_surface_correction(Dict("Surface Correction" => Dict("a" => 0)),
                                                                                                         "local_synch",
                                                                                                         "synchronise_field")
end
