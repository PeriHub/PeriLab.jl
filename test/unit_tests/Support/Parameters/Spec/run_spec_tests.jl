# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Standalone runner for the ParameterSpec unit tests:
#   julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl
using Test
import PeriLab

@testset "ParameterSpec" begin
    include(joinpath(@__DIR__, "spec_tests.jl"))
end
