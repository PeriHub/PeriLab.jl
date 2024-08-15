# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../helpers.jl")
using .Helpers: interpolation, interpol_data
export get_model_parameter
export find_data_files

"""
    get_model_parameter(params, model, id)

Retrieve a model parameter from a dictionary of parameters.

This function retrieves a specific model parameter from a dictionary of parameters based on the provided model and identifier (id).

## Arguments

- `params::Dict`: A dictionary containing various parameters.

- `model::String`: The model type for which the parameter is sought.

- `id::String`: The identifier (name) of the specific model parameter.

## Returns

- `parameter::Any`: The retrieved model parameter, or `nothing` if the parameter is not found.

## Errors

- If the specified model is defined in blocks but no model definition block exists, an error message is logged, and the function returns `nothing`.

- If the model with the given identifier is defined in blocks but missing in the model's definition, an error message is logged, and the function returns `nothing`.

## Example

```julia
params = Dict(
    "Models" => Dict(
        "Models" => Dict(
            "ModelA" => 42,
            "ModelB" => 24
        )
    )
)

model = "Models"
id = "ModelA"

result = get_model_parameter(params, model, id)
if result !== nothing
    println("Parameter id: result")
else
    println("Parameter not found.")
end
"""
function get_model_parameter(params::Dict, model::String, id::String)
    if !haskey(params["Models"], model * "s")
        @error model *
               " is defined in blocks, but no " *
               model *
               "s definition block exists"
        return nothing
    end
    if haskey(params["Models"][model*"s"], id)
        file_keys = find_data_files(params["Models"][model*"s"][id])
        for file_key in file_keys
            data = csv_reader_temporary(params["Models"][model*"s"][id][file_key])
            params["Models"][model*"s"][id][file_key] =
                interpolation(data[!, 1], data[!, 2])
        end
        return params["Models"][model*"s"][id]
    else
        @error model *
               " model with name " *
               id *
               " is defined in blocks, but missing in the " *
               model *
               "s definition."
        return nothing
    end
end
function csv_reader_temporary(filename::String)
    header_line, header = get_header(filename)
    return CSV.read(
        filename,
        DataFrame;
        delim = " ",
        ignorerepeated = true,
        header = header,
        skipto = header_line + 1,
        comment = "#",
    )
end


function find_data_files(params::Dict)
    file_keys = []
    for (key, value) in params
        if typeof(value) == String && endswith(value, ".txt")
            push!(file_keys, key)
        end
    end
    return file_keys
end


"""
    get_models_option(params, options)

Process models-related options based on the provided parameters.

This function processes models-related options based on the parameters dictionary and updates the options dictionary accordingly.

## Arguments

- `params::Dict`: A dictionary containing various parameters, including models-related information.

- `options::Dict`: A dictionary containing options to be updated based on the models parameters.

## Returns

- `updated_options::Dict`: A dictionary containing updated options based on the models parameters.

## Errors

- If the 'Pre Calculation' section exists in the 'Models' block but does not contain required options, an error message is logged.

- If a material model is missing the 'Material Model' specification, an error message is logged.

## Example

```julia
params = Dict(
    "Models" => Dict(
        "Pre Calculation" => Dict(
            "Option1" => true,
            "Option2" => false
        ),
        "Material Models" => Dict(
            1 => Dict(
                "Material Model" => "Correspondence"
            ),
            2 => Dict(
                "Material Model" => "Bond Associated"
            )
        )
    )
)

options = Dict(
    "Option1" => false,
    "Option2" => true
)

updated_options = get_models_option(params, options)
println("Updated Options: ", updated_options)
"""
function get_models_option(params::Dict, options::Dict)
    if haskey(params["Models"], "Pre Calculation")
        for option in keys(options)
            if haskey(params["Models"]["Pre Calculation"], option)
                options[option] = params["Models"]["Pre Calculation"][option]
            end
        end
    end
    if !haskey(params["Models"], "Material Models")
        @warn "Material Models missing!"
        return options
    end
    materials = params["Models"]["Material Models"]
    for material in eachindex(materials)
        if haskey(materials[material], "Material Model")
            options["Deformed Bond Geometry"] = true
            if occursin("Correspondence", materials[material]["Material Model"])
                if haskey(materials[material], "Bond Associated") &&
                   !(options["Bond Associated Deformation Gradient"])
                    # if its activated it stays that way
                    options["Bond Associated Deformation Gradient"] =
                        materials[material]["Bond Associated"]
                end
                if !(options["Bond Associated Deformation Gradient"])
                    options["Shape Tensor"] = true
                    options["Deformation Gradient"] = true
                end
            end
        else
            @error "No Material Model: '$material' has been defined"
            return nothing
        end

    end
    return options
end
