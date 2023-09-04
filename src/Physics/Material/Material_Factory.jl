
include("../../Support/tools.jl")
include("../../Support/Parameters/parameter_handling.jl")
module Material
export compute_forces

function compute_forces(datamanager, material, time, dt)
    # check which models
    if material["Material Model"] == "PD Solid Elastic"
        return Elastic.compute_forces(datamanager, material, time, dt)
    end
    @error "No material of type " * material["Material Model"] * " exists."
    return datamanager
end

function testing_material(params, datamanager)
    # for testing
    return datamanager
end
end