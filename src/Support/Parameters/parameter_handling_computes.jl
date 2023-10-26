# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_computes_names(params::Dict)
    if check_element(params::Dict, "Compute Class Parameters")
        computes = params["Compute Class Parameters"]
        return string.(collect(keys(sort(computes))))
    end
    return String[]
end

function get_output_variables(output, variables::Vector)
    if output in variables
        return output
    elseif output * "NP1" in variables
        return output * "NP1"
    else
        @warn '"' * output * '"' * " is not defined as variable"
    end
end

function get_computes(params::Dict, variables)
    if check_element(params::Dict, "Compute Class Parameters")
        computes = params["Compute Class Parameters"]
        for compute in keys(computes)
            if check_element(computes[compute], "Variable")
                computes[compute]["Variable"] = get_output_variables(computes[compute]["Variable"], variables)
            else
                @warn "No output variables are defined for " * output * "."
            end

        end
        return computes
    end
    return []
end

function get_node_set(params::Dict)
    if !check_element(params::Dict, "Node Set")
        return []
    end
    nodeset = params["Node Set"]

    if (typeof(nodeset) == Int64) | (typeof(nodeset) == Int32)
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