# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module Parameter_Handling
using Reexport
include("./parameter_handling_bc.jl")
include("./parameter_handling_blocks.jl")
include("./parameter_handling_physics.jl")
include("./parameter_handling_mesh.jl")
include("./parameter_handling_output.jl")
include("./parameter_handling_computes.jl")
include("./parameter_handling_solver.jl")
include("./parameter_handling_FEM.jl")

export validate_yaml

global expected_structure = Dict(
    "PeriLab" => [Dict{Any,Any}(
            "Blocks" => [Dict{Any,Any}(
                    "Any" => [Dict{Any,Any}(
                            "Density" => [Union{Float64,Int64}, true],
                            "Horizon" => [Union{Float64,Int64}, true],
                            "Specific Heat Capacity" => [Union{Float64,Int64}, false],
                            "Material Model" => [String, false],
                            "Damage Model" => [String, false],
                            "Thermal Model" => [String, false],
                            "Additive Model" => [String, false],
                        ), true],
                ), true],
            "FEM" => [Dict{Any,Any}(
                    "Element Type" => [String, true],
                    "Degree" => [Union{String,Int64}, true],
                    "Material Model" => [String, true],
                ), false],
            "Boundary Conditions" => [Dict{Any,Any}(
                    "Any" => [Dict{Any,Any}(
                            "Coordinate" => [String, false],
                            "Node Set" => [String, true],
                            "Type" => [String, true],
                            "Value" => [Union{Float64,Int64,String}, true],
                        ), true],
                ), false],
            "Compute Class Parameters" => [Dict{Any,Any}(
                    "Any" => [Dict{Any,Any}(
                            "Block" => [String, false],
                            "Node Set" => [String, false],
                            "Calculation Type" => [String, true],
                            "Compute Class" => [String, true],
                            "Variable" => [String, true],
                        ), true],
                ), false],
            "Discretization" => [Dict{Any,Any}(
                    "Input Mesh File" => [String, true],
                    "Input External Topology" => [Dict{Any,Any}(
                            "File" => [String, true],
                            "Add Neighbor Search" => [Bool, false],
                        ), false],
                    "Node Sets" => [Dict{Any,Any}(
                            "Any" => [Union{Int64,String}, true],
                        ), false],
                    "Type" => [String, true],
                    "Distribution Type" => [String, false],
                    "Surface Extrusion" => [Dict{Any,Any}(
                            "Direction" => [String, true],
                            "Step" => [Union{Float64,Int64}, true],
                            "Number" => [Union{Float64,Int64}, true],
                        ), false],
                    "Bond Filters" => [Dict{Any,Any}(
                            "Any" => [Dict{Any,Any}(
                                    "Type" => [String, true],
                                    "Normal X" => [Union{Float64,Int64}, true],
                                    "Normal Y" => [Union{Float64,Int64}, true],
                                    "Normal Z" => [Union{Float64,Int64}, true],
                                    "Lower Left Corner X" => [Union{Float64,Int64}, false],
                                    "Lower Left Corner Y" => [Union{Float64,Int64}, false],
                                    "Lower Left Corner Z" => [Union{Float64,Int64}, false],
                                    "Bottom Unit Vector X" => [Union{Float64,Int64}, false],
                                    "Bottom Unit Vector Y" => [Union{Float64,Int64}, false],
                                    "Bottom Unit Vector Z" => [Union{Float64,Int64}, false],
                                    "Center X" => [Union{Float64,Int64}, false],
                                    "Center Y" => [Union{Float64,Int64}, false],
                                    "Center Z" => [Union{Float64,Int64}, false],
                                    "Radius" => [Union{Float64,Int64}, false],
                                    "Bottom Length" => [Union{Float64,Int64}, false],
                                    "Side Length" => [Union{Float64,Int64}, false],
                                    "Allow Contact" => [Bool, false],
                                ), true],
                        ), false],
                ), true],
            "Outputs" => [Dict{Any,Any}(
                    "Any" => [Dict{Any,Any}(
                            "Flush File" => [Bool, false],
                            "Output Frequency" => [Int64, false],
                            "Number of Output Steps" => [Int64, false],
                            "Output File Type" => [String, false],
                            "Output Filename" => [String, true],
                            "Write After Damage" => [Bool, false],
                            "Output Variables" => [Dict{Any,Any}(
                                    "Any" => [Bool, true],
                                ), true],
                        ), true],
                ), false],
            "Physics" => [Dict{Any,Any}(
                    "Damage Models" => [Dict{Any,Any}(
                            "Any" => [Dict{Any,Any}(
                                    "Critical Value" => [Union{Float64,Int64}, true],
                                    "Damage Model" => [String, true],
                                    "Interblock Damage" => [Dict{Any,Any}(
                                            "Any" => [Union{Float64,Int64}, true],
                                        ), false],
                                    "Anisotropic Damage" => [Dict{Any,Any}(
                                            "Critical Value X" => [Union{Float64,Int64}, true],
                                            "Critical Value Y" => [Union{Float64,Int64}, true],
                                        ), false],
                                ), true],
                        ), false],
                    "Material Models" => [Dict{Any,Any}(
                            "Any" => [Dict{Any,Any}(
                                    "Material Model" => [String, true],
                                    "Symmetry" => [String, false],
                                    "Poisson's Ratio" => [Union{Float64,Int64}, false],
                                    "Young's Modulus" => [Union{Float64,Int64}, false],
                                    "Bulk Modulus" => [Union{Float64,Int64}, false],
                                    "Shear Modulus" => [Union{Float64,Int64}, false],
                                    "Yield Stress" => [Union{Float64,Int64}, false],
                                    "Zero Energy Control" => [String, false],
                                    "C11" => [Union{Float64,Int64}, false],
                                    "C12" => [Union{Float64,Int64}, false],
                                    "C13" => [Union{Float64,Int64}, false],
                                    "C14" => [Union{Float64,Int64}, false],
                                    "C15" => [Union{Float64,Int64}, false],
                                    "C16" => [Union{Float64,Int64}, false],
                                    "C22" => [Union{Float64,Int64}, false],
                                    "C23" => [Union{Float64,Int64}, false],
                                    "C24" => [Union{Float64,Int64}, false],
                                    "C25" => [Union{Float64,Int64}, false],
                                    "C26" => [Union{Float64,Int64}, false],
                                    "C33" => [Union{Float64,Int64}, false],
                                    "C34" => [Union{Float64,Int64}, false],
                                    "C35" => [Union{Float64,Int64}, false],
                                    "C36" => [Union{Float64,Int64}, false],
                                    "C44" => [Union{Float64,Int64}, false],
                                    "C45" => [Union{Float64,Int64}, false],
                                    "C46" => [Union{Float64,Int64}, false],
                                    "C55" => [Union{Float64,Int64}, false],
                                    "C56" => [Union{Float64,Int64}, false],
                                    "C66" => [Union{Float64,Int64}, false],
                                ), true],
                        ), false],
                    "Thermal Models" => [Dict{Any,Any}(
                            "Any" => [Dict{Any,Any}(
                                    "Thermal Model" => [String, true],
                                    "Type" => [String, false],
                                    "Heat Transfer Coefficient" => [Union{Float64,Int64}, false],
                                    "Environmental Temperature" => [Union{Float64,Int64}, false],
                                    "Thermal Conductivity" => [Union{Float64,Int64}, false],
                                    "Thermal Expansion Coefficient" => [Union{Float64,Int64}, false],
                                    "Thermal Conductivity Print Bed" => [Union{Float64,Int64}, false],
                                    "Print Bed Temperature" => [Union{Float64,Int64}, false],
                                ), true],
                        ), false],
                    "Additive Models" => [Dict{Any,Any}(
                            "Any" => [Dict{Any,Any}(
                                    "Additive Model" => [String, true],
                                    "Print Temperature" => [Union{Float64,Int64}, false],
                                ), true],
                        ), false],
                    "Pre Calculation" => [Dict{Any,Any}(
                            "Bond Associated Deformation Gradient" => [Bool, false],
                            "Bond Associated Shape Tensor" => [Bool, false],
                            "Deformation Gradient" => [Bool, false],
                            "Deformed Bond Geometry" => [Bool, false],
                            "Shape Tensor" => [Bool, false],
                        ), false],
                ), true],
            "Solver" => [Dict{Any,Any}(
                     "Material Models" => [Bool, false],
                    "Damage Models" => [Bool, false],
                    "Maximum Damage" => [Float64, false],
                    "Thermal Models" => [Bool, false],
                    "Additive Models" => [Bool, false],
                    "Final Time" => [Union{Float64,Int64}, true],
                    "Initial Time" => [Union{Float64,Int64}, true],
                    "Numerical Damping" => [Union{Float64,Int64}, false],
                    "Verlet" => [Dict{Any,Any}(
                            "Safety Factor" => [Union{Float64,Int64}, false],
                            "Fixed dt" => [Union{Float64,Int64}, false],
                            "Number of Steps" => [Int64, false]
                        ), false],
                    "External" => [Dict{Any,Any}(
                            "Number of Steps" => [Int64, false]
                        ), false]
                ), true],
        ), true]
)

"""
    validate_structure_recursive(expected::Dict, actual::Dict, validate::Bool, checked_keys::Array, path::String="")

Validates the parameters against the expected structure

# Arguments
- `expected::Dict`: The expected structure
- `actual::Dict`: The actual structure
- `validate::Bool`: The validation results
- `checked_keys::Array`: The keys that have been checked
- `path::String`: The current path
# Returns
- `validate::Bool`: The validation result
- `checked_keys::Array`: The keys that have been checked
"""
function validate_structure_recursive(expected::Dict, actual::Dict, validate::Bool, checked_keys::Array, path::String="")
    for (key, value) in expected
        current_path = isempty(path) ? key : "$path.$key"

        if key == "Any"
            if length(keys(actual)) == 0 && value[2] == true
                @error "Validation Error: Missing key - $current_path"
                validate = false
                continue
            end
            for (any_key, any_value) in actual
                if isa(any_value, Dict) && isa(actual[any_key], Dict)
                    # Recursive call for nested dictionaries
                    validate, checked_keys = validate_structure_recursive(expected[key][1], actual[any_key], validate, checked_keys, current_path)
                end
                push!(checked_keys, any_key)
            end
            continue
        end

        if !haskey(actual, key) && value[2] == true
            @error "Validation Error: Missing key - $current_path"
            validate = false
            continue
        end

        if !haskey(actual, key) && value[2] == false
            continue
        end

        if isa(actual[key], typeof(value[1])) || isa(actual[key], value[1])
            push!(checked_keys, key)
            if isa(value[1], Dict) && isa(actual[key], Dict)
                # Recursive call for nested dictionaries
                validate, checked_keys = validate_structure_recursive(value[1], actual[key], validate, checked_keys, current_path)
            end
        else
            @error "Validation Error: Wrong type, expected - $(value[1]), got - $(typeof(actual[key])) in $current_path"
            validate = false
        end
    end
    return validate, checked_keys
end

"""
    get_all_keys(params::Dict)

Get all the keys in the parameters

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `keys_list::Array`: The keys list
"""
function get_all_keys(params::Dict)
    keys_list = []
    for (key, value) in params
        push!(keys_list, key)
        if isa(value, Dict)
            keys_list = vcat(keys_list, get_all_keys(value))
        end
    end
    return keys_list
end

"""
    validate_yaml(params::Dict)

Validates the parameters against the expected structure

# Arguments
- `params::Dict`: The parameters dictionary.
# Returns
- `params::Dict`: The parameters dictionary.
"""
function validate_yaml(params::Dict)
    all_keys = get_all_keys(params)
    # Validate against the expected structure
    validate = true
    checked_keys = []
    validate, checked_keys = validate_structure_recursive(expected_structure, params, validate, checked_keys)
    #Check if all keys have been checked
    for key in all_keys
        if !(key in checked_keys)
            @warn "Key not known - $key, going to ignore it"
        end
    end
    if !validate
        return nothing
    end

    return params["PeriLab"]
end
end