
include("../../Support/tools.jl")
include("../../Support/Parameters/parameter_handling.jl")
module Material
include("./Ordinary/Ordinary.jl")
using .PD_Solid_Elastic
export compute_forces

function compute_forces(datamanager, material, time, dt)
    # check which models
    if material["Material Model"] == "Test"
        return testing_material(datamanager, time)
    end

    if material["Material Model"] == "PD Solid Elastic"
        return PD_Solid_Elastic.compute_forces(datamanager, material, time, dt)
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