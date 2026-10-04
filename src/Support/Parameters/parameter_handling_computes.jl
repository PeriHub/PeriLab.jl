# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using DataStructures: OrderedDict

export get_computes_names
export get_output_variables
export get_computes
export compute_names, active_computes

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
            output = get_output_variables(params["Compute Class Parameters"][compute]["Variable"],
                                          variables)
            if isnothing(output)
                continue
            end
            computes[compute] = params["Compute Class Parameters"][compute]
            computes[compute]["Variable"] = output
        else
            @warn "No output variables are defined for " *
                  output *
                  ". Global variable is not defined"
        end
    end
    return computes
end

# Typed functions (phase 2b). The Dict getters above stay until phase 4.

"Names of the compute classes, sorted."
compute_names(computes::Dict{String,ComputeClassParams}) = sort!(collect(keys(computes)))

"""
    active_computes(computes, field_keys)

The compute classes whose `Variable` is an existing field (directly or as
`<Variable>NP1`); the others are skipped with a warning.
"""
function active_computes(computes::Dict{String,ComputeClassParams},
                         field_keys::Vector{String})
    active = Dict{String,ComputeClassParams}()
    for (name, compute) in computes
        if compute.variable in field_keys || compute.variable * "NP1" in field_keys
            active[name] = compute
        else
            @warn '"' * compute.variable * '"' * " is not defined as variable"
        end
    end
    return active
end
