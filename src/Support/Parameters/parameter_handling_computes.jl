# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_computes_names(params::Dict)
    if check_element(params::Dict, "Compute Class Parameters")
        computes = params["Compute Class Parameters"]
        return collect(keys(sort(computes)))
    end
    return []
end

function get_output_variables(output, variables)
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