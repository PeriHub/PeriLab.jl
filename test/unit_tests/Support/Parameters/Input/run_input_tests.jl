# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Standalone runner for the InputDeck unit tests:
#   julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl
using Test
import PeriLab
include(joinpath(@__DIR__, "..", "..", "..", "..", "helper.jl"))

@testset "InputDeck" begin
    include(joinpath(@__DIR__, "input_tests.jl"))
end
