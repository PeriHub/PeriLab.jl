# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

for file in ["ut_errors.jl", "ut_dependent.jl", "ut_convert.jl", "ut_params.jl", "ut_model.jl",
             "ut_end_to_end.jl", "ut_extensions.jl"]
    @testset "$file" begin
        include(joinpath(@__DIR__, file))
    end
end
