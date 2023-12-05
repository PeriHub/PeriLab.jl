module FEM
include("../Core/Module_inclusion/set_Modules.jl")
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


using .Set_modules

#global module_list = Set_modules.find_module_files(@__DIR__, "element_name")
#Set_modules.include_files(module_list)

#datamanager = Set_modules.create_module_specifics(fem_model, module_list, specifics, (datamanager, nodes, model_param, time, dt))
#if isnothing(datamanager)
#    @error "No shape function of name " * fem_model * " exists."
#end


function get_polynomial_degree(params::Dict, dof::Int64)
    if !haskey(params, "Degree")
        @error "No element degree defined"
        return nothing
    end
    value = params["Degree"]
    if sum(typeof.(value) .!= Int64) != 0
        @warn "Degree was defined as Float and set to Int."
        value = Int64.(round.(value))
    end
    if length(value) == 1
        return_value::Vector{Int64} = zeros(dof)
        return_value[1:dof] .= value[1]
        return return_value
    elseif length(value) == dof
        return value[1:dof]
    else
        @error "Degree must be defined with length one or number of dof."
        return nothing
    end
end


end