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

    # no section: nothing is stored and the compute step is a no-op
    @test isnothing(mod_struct.init_surface_correction(nothing, "local_synch",
                                                       "synchronise_field"))
    @test isnothing(PeriLab.Data_Manager.get_surface_correction())
    @test isnothing(mod_struct.compute_surface_correction([1, 2], "local_synch",
                                                          "synchronise_field"))

    sc = typed_section(PeriLab.InputDeck.SurfaceCorrectionParams,
                       Dict{String,Any}("Type" => "Volume Correction", "Update" => true))
    PeriLab.Data_Manager.set_surface_correction(sc)
    @test PeriLab.Data_Manager.get_surface_correction().update
end
