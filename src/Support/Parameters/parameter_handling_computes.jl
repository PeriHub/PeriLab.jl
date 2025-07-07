# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using OrderedCollections: OrderedDict

export get_computes_names
export get_output_variables
export get_computes

"""
    get_computes_names(params::Dict)

Get the names of the computes.

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `computes_names::Vector{String}`: The names of the computes.
"""
function get_computes_names(params::Dict)
    if haskey(params::Dict, "Compute Class Parameters")
        computes = params["Compute Class Parameters"]
        return string.(collect(keys(sort!(OrderedDict(computes)))))
    end
    return String[]
end

"""
    get_output_variables(output::String, variables::Vector{String})

Get the output variable.

# Arguments
- `output::String`: The output variable.
- `variables::Vector{String}`: The variables.
# Returns
- `output::String`: The output variable.
"""
function get_output_variables(output::String, variables::Vector{String})
    if output in variables || output * "NP1" in variables
        return output
    else
        @warn '"' * output * '"' * " is not defined as variable"
    end
end

"""
    get_computes(params::Dict, variables::Vector{String})

Get the computes.

# Arguments
- `params::Dict`: The parameters dictionary.
- `variables::Vector{String}`: The variables.
# Returns
- `computes::Dict{String,Dict{Any,Any}}`: The computes.
"""
function get_computes(params::Dict, variables::Vector{String})
    computes = Dict{String,Dict{Any,Any}}()
    if !haskey(params, "Compute Class Parameters")
        return computes
    end
    for compute in keys(params["Compute Class Parameters"])
        if haskey(params["Compute Class Parameters"][compute], "Variable")
            computes[compute] = params["Compute Class Parameters"][compute]
            computes[compute]["Variable"] = get_output_variables(computes[compute]["Variable"],
                                                                 variables)
        else
            @warn "No output variables are defined for " *
                  output *
                  ". Global variable is not defined"
        end
    end
    return computes
end
