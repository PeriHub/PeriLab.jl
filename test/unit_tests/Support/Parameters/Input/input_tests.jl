# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl",
             "ut_read_input.jl", "ut_golden_decks.jl", "ut_validate_yaml.jl",
             "ut_solver_typed.jl", "ut_model_reduction_params.jl", "ut_blocks_mesh_typed.jl",
             "ut_node_sets_typed.jl", "ut_outputs_typed.jl"]
    @testset "$file" begin
        include(joinpath(@__DIR__, file))
    end
end
