# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_model_parameter(params::Dict, model::String, id::Int64)
    if check_element(params["Physics"], model * "s") == false
        @error model * " is defined in blocks, but no " * model * "s definition block exists"
        return Dict()
    end
    if check_element(params["Physics"][model*"s"], id)
        return params["Physics"][model*"s"][id]
    else
        @error model * " model with name " * id * " is defined in blocks, but missing in the " * model * "s definition."
        return Dict()
    end
end

function get_physics_option(params::Dict, options::Dict)
    if check_element(params["Physics"], "Pre Calculation")
        for option in keys(options)
            if check_element(params["Physics"]["Pre Calculation"], option)
                options[option] = params["Physics"]["Pre Calculation"][option]
            end
        end
    end
    if !haskey(params["Physics"], "Material Models")
        @error "Material Models missing!"
    end
    materials = params["Physics"]["Material Models"]
    for material in eachindex(materials)
        if haskey(materials[material], "Material Model")
            if occursin("Correspondence", materials[material]["Material Model"])
                options["Shape Tensor"] = true
                options["Deformation Gradient"] = true
                options["Deformed Bond Geometry"] = true
            elseif occursin("Bond Associated", materials[material]["Material Model"])
                options["Shape Tensor"] = true
                options["Bond Associated Shape Tensor"] = true
                options["Bond Associated Deformation Gradient"] = true
            end
        else
            @error "Material $material is missing 'Material Model'!"
        end
    end
    return options
end

