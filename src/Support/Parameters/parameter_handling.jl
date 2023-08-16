include("./parameter_handling_bc.jl")
include("./parameter_handling_blocks.jl")
include("./parameter_handling_physics.jl")
include("./parameter_handling_mesh.jl")
include("./parameter_handling_output.jl")
include("./parameter_handling_solver.jl")

function check_element(params, key)
    return haskey(params, key)
end
function check_key_elements(params)
    check = check_element(params, "Material Models")
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
