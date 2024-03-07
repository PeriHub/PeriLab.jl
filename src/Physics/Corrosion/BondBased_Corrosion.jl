# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Corrision
export compute_corrosion_model
export corrosion_name
export init_corrosion_model

"""
    corrosion_name()

Gives the corrosion name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the bond-based corrosion model.

Example:
```julia
println(corrosion_name())
"Bond-based Corrosion"
```
"""
function corrosion_name()
    return "Bond-based Corrosion"
end

"""
    compute_corrosion_model(datamanager, nodes, corrosion_parameter, time, dt)

Calculates the bond-based corrosion model. This template has to be copied, the file renamed and edited by the user to create a new corrosion. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `corrosion parameter::Dict(String, Any)`: Dictionary with corrosion parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
  ```
"""
function compute_corrosion_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, corrosion_parameter::Dict, time::Float64, dt::Float64)
    concentrationN = datamanager.get_field("Concentration", "N")
    concentrationNP1 = datamanager.get_field("Concentration", "NP1")

    return datamanager
end


"""
    init_corrosion_model(datamanager, nodes, corrosion_parameter)

Inits the bond-based corrosion model. This template has to be copied, the file renamed and edited by the user to create a new corrosion. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `corrosion parameter::Dict(String, Any)`: Dictionary with corrosion parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_corrosion_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, corrosion_parameter::Dict, block::Int64)
 
    return datamanager
end


end