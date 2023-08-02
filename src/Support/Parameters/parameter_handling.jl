include("./parameter_handling_bc.jl")
include("./parameter_handling_blocks.jl")
include("./parameter_handling_material.jl")
include("./parameter_handling_mesh.jl")
include("./parameter_handling_output.jl")
include("./parameter_handling_solver.jl")

function check_element(params, key)
    return haskey(params, key)
end
function check_key_elements(params)
    check = check_element(params, "Materials")
    if !check
        @error "No materials defined"
        return check
    end
    check = check_element(params, "Blocks")
    if !check
        @error "No blocks defined"
        return check
    end
    check = check_element(params, "Discretization")
    if !check
        @error "No discretization defined"
        return check
    end
    check = check_element(params, "Solver")
    if !check
        @error "No solver defined"
        return check
    end
    return check
end
function get_element(params, key)
    check = check_element(params, "Materials")
    if !check
        @error "Entry " * key * " is not defined"
        return check
    end
    return params[key]
end