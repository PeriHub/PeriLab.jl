# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    ParameterSpec

Typed, validated input parameters. Modules declare their parameters with
`@params`; YAML input is converted into those structs and validated against
the declarations. See
`docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md`.
"""
module ParameterSpec

using ..PeriLabExceptions: @abort
using Dierckx: Spline1D, evaluate

include("errors.jl")
include("dependent.jl")
include("field_spec.jl")
include("convert.jl")
include("params_macro.jl")
include("build.jl")
include("registry.jl")
include("model.jl")
include("bind.jl")

end
