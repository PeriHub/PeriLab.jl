# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Corrosion

using .......Data_Manager
export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_degradation
@params struct BondbasedCorrosionParams
end
__init__() = register_degradation("Bond-based Corrosion", BondbasedCorrosionParams)

"""
    compute_model(nodes, p, degradation, block, time, dt)

Calculates the bond-based degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::BondbasedCorrosionParams`: The model parameters.
- `degradation`: `nothing` (degradation models have no shared parameters).
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
  ```
"""
function compute_model(nodes::AbstractVector{Int64}, p::BondbasedCorrosionParams,
                       degradation, block::Int64, time::Float64, dt::Float64)
    concentrationN = Data_Manager.get_field("Concentration", "N")
    concentrationNP1 = Data_Manager.get_field("Concentration", "NP1")
    concentration_fluxN = Data_Manager.get_field("Concentration Flux", "N")
    concentration_fluxNP1 = Data_Manager.get_field("Concentration Flux", "NP1")
end

"""
    init_model(nodes, p, degradation, block)

Inits the bond-based degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::BondbasedCorrosionParams`: The model parameters.
- `degradation`: `nothing` (degradation models have no shared parameters).
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::BondbasedCorrosionParams, degradation,
                    block::Int64)
end

"""
    fields_for_local_synchronization(dmodel::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
