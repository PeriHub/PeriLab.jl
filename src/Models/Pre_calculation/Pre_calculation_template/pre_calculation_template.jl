# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_calculation_template

using .......Data_Manager
export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_pre_calculation
"Switch only: this pre-calculation has no parameters."
@params struct PreCalculationTemplateParams
end
# __init__() = register_pre_calculation("pre_calculation Template", PreCalculationTemplateParams)
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
    compute_model(nodes, p, pre_calculation, block, time, dt)

This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::PreCalculationTemplateParams`: The model parameters.
- `pre_calculation`: `nothing` (pre-calculations have no shared parameters).
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
  ```
"""
function compute_model(nodes::AbstractVector{Int64}, p::PreCalculationTemplateParams,
                       pre_calculation, block::Int64, time::Float64, dt::Float64)
    @info "Register your model name with register_pre_calculation in __init__()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(nodes, p, pre_calculation, block, time, dt) function."
    @info "The Data_Manager holds all you need to solve your problem on material level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

"""
    init_model(nodes, p, pre_calculation, block)

Inits the calculation.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::PreCalculationTemplateParams`: The model parameters.
- `pre_calculation`: `nothing` (pre-calculations have no shared parameters).
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::PreCalculationTemplateParams,
                    pre_calculation, block::Int64)
end

end
