# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_calculation_template

using .......Data_Manager
export compute
export init_model
export pre_calculation_name
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_pre_calculation
"Switch only: this pre-calculation has no parameters."
@params struct PreCalculationTemplateParams
end
# __init__() = register_pre_calculation("pre_calculation Template", PreCalculationTemplateParams)
"""
    pre_calculation_name()

Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the Pre_Calculation.

Example:
```julia
println(pre_calculation_name())
"Pre_calculation Template"
```
"""
function pre_calculation_name()
    return "pre_calculation Template"
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

"""
    compute(nodes, block)

This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `block::Int64`: The current block.
Example:
```julia
  ```
"""
function compute(nodes::AbstractVector{Int64}, block::Int64)
    @info "Please write a possible precalculation routines in pre_calculation_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute(nodes, block) function."
    @info "The Data_Manager holds all you need to solve your problem on material level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

"""
    init_model(nodes, block)

Inits the calculation.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.

"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
end

end
