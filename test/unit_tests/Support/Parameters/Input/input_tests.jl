# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl"]
    @testset "$file" begin
        include(joinpath(@__DIR__, file))
    end
end
