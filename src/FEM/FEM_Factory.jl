module FEM
include("../Core/Module_inclusion/set_Modules.jl")
include("./FEM_routines.jl")
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "element_name")
Set_modules.include_files(module_list)

function init_FEM(datamanager::Module, params::Dict)
    valid_models(params)
    dof = datamanager.get_dof()
    nelements = datamanager.get_num_elements()
    elements::Vector{Int64} = 1:nelements
    p = get_polynomial_degree(params, dof)

    if isnothing(p)
        return p
    end
    if dof != 2 && dof != 3
        @error "Degree of freedom = $dof is not supported, only 2 and 3."
        return nothing
    end
    num_int = get_number_of_integration_points(p, dof)
    N = datamanager.create_constant_free_size_field("N Matrix", Float64, (prod(num_int), prod(p .+ 1) * dof, dof))
    B = datamanager.create_constant_free_size_field("B Matrix", Float64, (prod(num_int), prod(p .+ 1) * dof, 3 * dof - 3))

    strainN, strainNP1 = datamanager.create_free_size_field("Element Strain", Float64, (nelements, prod(num_int), 3 * dof - 3))
    stressN, stressNP1 = datamanager.create_free_size_field("Element Stress", Float64, (nelements, prod(num_int), 3 * dof - 3))
    strain_increment = datamanager.create_constant_free_size_field("Element Strain Increment", Float64, (nelements, prod(num_int), 3 * dof - 3))
    if isnothing(N) || isnothing(B)
        return nothing
    end
    specifics = Dict{String,String}("Call Function" => "create_element_matrices", "Name" => "element_name")
    N[:], B[:] = create_element_matrices(dof, p, Set_modules.create_module_specifics(params["Element Type"], module_list, specifics))

    specifics = Dict{String,String}("Call Function" => "init_element", "Name" => "element_name")
    datamanager = Set_modules.create_module_specifics(params["Element Type"], module_list, specifics, (datamanager, elements, params, p))

    return datamanager

end



function valid_models(params::Dict)
    if haskey(params, "Additive Model")
        @warn "Additive models are not supported for FEM yet"
    end
    if haskey(params, "Damage Model")
        @warn "Damage models are not supported for FEM"

    end
    if haskey(params, "Thermal Model")
        @warn "Thermal models are not supported for FEM yet"
    end
    if !haskey(params, "Material Model")
        @error "No material model has been defined for the block"
        return nothing
    else
        if !occursin(params["Material Model"], "Correspondence")
            @error "Only correspondence material is supported for FEM"
            return nothing
        end
    end
    return nothing
end
function eval(datamanager::Module, elements::Union{SubArray,Vector{Int64}})
    return calculate_FEM(datamanager, elements)
end

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