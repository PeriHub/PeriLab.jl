include("./Parameters/parameter_mesh_handling.jl")
include("./Parameters/parameter_handling_blocks.jl")
function check_element(data, key)
    return haskey(data, key)
end