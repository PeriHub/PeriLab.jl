# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_model_parameter(params, model, id)
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

function get_physics_options(params, options)
    if check_element(params["Physics"], "Pre Calculation")
        for option in keys(options)
            if check_element(params["Physics"]["Pre Calculation"], option)
                options[option] = params["Physics"]["Pre Calculation"][option]
            end
        end
    end
    materials = params["Physics"]["Material Models"]
    for material in eachindex(materials)
        if occursin("Correspondence", materials[material]["Material Model"])
            options["Shape Tensor"] = true
            options["Deformation Gradient"] = true
            options["Deformed Bond Geometry"] = true
        end
        if occursin("Bond Associated", materials[material]["Material Model"])
            options["Shape Tensor"] = true
            options["Bond Associated Shape Tensor"] = true
            options["Bond Associated Deformation Gradient"] = true
        end
    end
    return options
end

