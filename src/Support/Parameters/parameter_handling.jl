# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("./parameter_handling_bc.jl")
include("./parameter_handling_blocks.jl")
include("./parameter_handling_physics.jl")
include("./parameter_handling_mesh.jl")
include("./parameter_handling_output.jl")
include("./parameter_handling_computes.jl")
include("./parameter_handling_solver.jl")

function check_element(params::Dict, key)
    return haskey(params::Dict, key)
end
function check_key_elements(params::Dict)
    if !check_element(params, "Physics")
        @error "No physics defined"
        return nothing
    end
    if length(params["Physics"]) == 0
        @error "No physics defined"
        return nothing
    end

    if !check_element(params, "Blocks")
        @error "No blocks defined"
        return nothing
    end

    if !check_element(params, "Discretization")
        @error "No discretization defined"
        return nothing
    end

    if !check_element(params, "Solver")
        @error "No solver defined"
        return nothing
    end
    return params
end
