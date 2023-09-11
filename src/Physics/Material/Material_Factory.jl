
include("../../Support/tools.jl")
include("../../Support/Parameters/parameter_handling.jl")
module Material
include("./Ordinary/PD_Solid_Elastic.jl")
using .PD_Solid_Elastic
export compute_forces

function compute_forces(datamanager, nodes, material, time, dt)
    # check which models
    if length(material) == 0
        @error "No material of type " * material["Material Model"] * " exists."
        return datamanager
    end
    if material["Material Model"] == "Test"
        return testing_material(datamanager, time)
    end

    if material["Material Model"] == PD_Solid_Elastic.material_name()
        return PD_Solid_Elastic.compute_forces(datamanager, nodes, material, time, dt)
    end
    @error "No material of type " * material["Material Model"] * " exists."
    return datamanager
end

function testing_material(datamanager, time)
    if time == 0
        @info "Testing dummy material is used. Nothing happens here"
    end
    # for testing
    return datamanager
end
end