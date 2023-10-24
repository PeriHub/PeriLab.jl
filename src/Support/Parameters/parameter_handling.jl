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
    check = check_element(params::Dict, "Physics")
    if !check
        @error "No physics defined"
        return check
    end
    if check
        if length(params["Physics"]) == 0
            @error "No physics defined"
            return false
        end
    end
    check = check_element(params::Dict, "Blocks")
    if !check
        @error "No blocks defined"
        return check
    end
    check = check_element(params::Dict, "Discretization")
    if !check
        @error "No discretization defined"
        return check
    end
    check = check_element(params::Dict, "Solver")
    if !check
        @error "No solver defined"
        return check
    end
    return check
end
