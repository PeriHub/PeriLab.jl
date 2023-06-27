include("./Parameters/parameter_mesh_handling.jl")

function check_element(data, key)
    return haskey(data, key)
end