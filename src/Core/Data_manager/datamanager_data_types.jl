# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# =============================================================================
# ABSTRACT BASE TYPES
# =============================================================================

abstract type DataField{T} end

include("./Types/NodeField.jl")
include("./Types/BondField.jl")
include("./Types/FreeSizeField.jl")
include("./Types/ElementField.jl")
