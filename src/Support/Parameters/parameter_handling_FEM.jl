# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../helpers.jl")


function get_element_type(params::Dict)
    if !haskey(params, "Element Type")
        @error "Element Type is not defined."
        return nothing
    end
    if typeof(params["Element Type"]) != String
        @warn "Element Type must be of type string."
    end
    return string(params["Element Type"])
end