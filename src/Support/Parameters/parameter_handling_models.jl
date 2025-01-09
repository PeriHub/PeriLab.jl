# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../Helpers.jl")
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
function get_model_parameter(
    params::Dict,
    model::String,
    id::String,
    directory::String = "",
)
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
            data, header = csv_reader_temporary(
                joinpath(directory, params["Models"][model*"s"][id][file_key]),
            )
            params["Models"][model*"s"][id][file_key] = Dict()
            params["Models"][model*"s"][id][file_key]["Field"] = header[1]
            params["Models"][model*"s"][id][file_key]["Data"] =
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
    ),
    header
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
