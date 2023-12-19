# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

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
        return string.(collect(keys(sort(computes))))
    end
    return String[]
end

""" 
    get_output_variables(output::String, variables::Vector)

Get the output variable.

# Arguments
- `output::String`: The output variable.
- `variables::Vector`: The variables.
# Returns
- `output::String`: The output variable.
"""
function get_output_variables(output::String, variables::Vector)
    if output in variables
        return output
    elseif output * "NP1" in variables
        return output * "NP1"
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
            computes[compute]["Variable"] = get_output_variables(computes[compute]["Variable"], variables)
        else
            @warn "No output variables are defined for " * output * ". Global variable is not defined"
        end
    end
    return computes
end

"""
    get_node_set(computes::Dict, path::String, params::Dict)

Get the node set.

# Arguments
- `computes::Dict`: The computes dictionary.
- `path::String`: The path to the mesh.
- `params::Dict`: The parameters dictionary.
# Returns
- `nodeset::Vector`: The node set.
"""
function get_node_set(computes::Dict, path::String, params::Dict)
    if !haskey(computes::Dict, "Node Set")
        return []
    end
    nodeset = computes["Node Set"]

    if params["Discretization"]["Type"] == "Exodus"
        exo = ExodusDatabase(joinpath(path, get_mesh_name(params)), "r")
        nset = read_set(exo, NodeSet, nodeset)
        return Vector{Int64}(nset.nodes)
    end

    if nodeset isa Int64 || nodeset isa Int32
        return [nodeset]
    elseif occursin(".txt", nodeset)

        header_line, header = get_header(nodeset)
        nodes = CSV.read(nodeset, DataFrame; delim=" ", header=false, skipto=header_line + 1)
        if size(nodes) == (0, 0)
            @error "Node set file is empty " * nodeset
        end
        return nodes.Column1
    else
        nodes = split(nodeset)
        return parse.(Int, nodes)
    end
end