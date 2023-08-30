
include("../../Support/tools.jl")
include("../../Support/Parameters/parameter_handling.jl")
module Material
export get_material


function material_types(params, datamanager)
    material_types = get_material_types(params)
    for material_type in material_types
        datamanager.set_material_type(material_type, material_types[material_type])
    end
    return datamanager.get_material_type()
end

function get_material(params)
    materials = get_materials(params)
    material_names = keys(material)
    #files = find_files_with_ending("./", "*.jl")
    # blocks
    # blocks -> mat
end

function evaluate_material(params, datamanager, dt, time)

    #datamanager.get
    #    if 
end


function create_nodel_forces()
    println()
end

function testing_material(params, datamanager)
    # for testing
    return datamanager
end
end