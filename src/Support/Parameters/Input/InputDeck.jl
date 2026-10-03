# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    InputDeck

Typed declarations of PeriLab's fixed input sections and `read_input`, which
turns the `PeriLab:` part of a YAML deck into a validated `PeriLabInput`.
Model parameters (`Models`) are declared by the model modules (phase 3).
"""
module InputDeck

using ..ParameterSpec
using ..ParameterSpec: @params, req, opt, ParseContext, add_error!, join_path,
                       parse_section
import ..ParameterSpec: check!

include("discretization.jl")
include("blocks.jl")
include("solver.jl")
include("outputs.jl")
include("conditions.jl")
include("contact.jl")
include("input.jl")

export read_input, PeriLabInput

end
