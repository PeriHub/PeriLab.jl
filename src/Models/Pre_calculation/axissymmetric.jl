# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Axissymmetric
export compute
export Pre_calculation_name

"""
    pre_calculation_name()

Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the Pre_calculation.

Example:
```julia
println(pre_calculation_name())
"Axis Symmetric"
```
"""
function pre_calculation_name()
    return "Axis Symmetric"
end

"""
    init_model(datamanager, nodes, corrosion_parameter)

Inits the bond-based corrosion model. This template has to be copied, the file renamed and edited by the user to create a new corrosion. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `corrosion parameter::Dict(String, Any)`: Dictionary with corrosion parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    corrosion_parameter::Dict,
    block::Int64,
)

    return datamanager
end


"""
    compute(datamanager, nodes, Pre_calculation_parameter, time, dt)

This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `Pre_calculation_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
  ```
"""
function compute(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    Pre_calculation_parameter::Dict,
    time::Float64,
    dt::Float64,
)
    @info "Please write a possible precalculation routines in pre_calculation_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(datamanager, nodes, Pre_calculation_parameter, time, dt) function."
    @info "The datamanager and Pre_calculation_parameter holds all you need to solve your problem on material level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end

function init(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    Pre_calculation_parameter::Dict,
)
    symmetry_axis = datamanager.get_symmetry_axis()
    volume = datamanager.get_field("Volume")
    coordinates = datamanager.get_field("Coordinates")
    # Volume must be an area
    for iID in nodes
        if coordinate[iID, 2] - symmetry_axis == 0
            volume[iID] *= 0.5
            continue
        end
        volume[iID] *= 2 * pi * coordinates[iID, 2] - symmetry_axis
    end
    return datamanager
end
end
