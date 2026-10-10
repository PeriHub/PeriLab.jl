# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Rigid

using .......Data_Manager
using ......ParameterSpec: @params, register_material

export init_model
export fe_support
export compute_model

"Parameters of Rigid beyond the shared material keys (none)."
@params struct RigidParams
end
__init__() = register_material("Rigid", RigidParams)

"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
    return true
end

"""
  init_model(nodes::AbstractVector{Int64}, p, material, block::Int64)

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `p`: Model parameters; `material::BlockMaterial`: typed block material (`material.base`, `material.moduli`, `material.symmetry`).
"""
function init_model(nodes::AbstractVector{Int64},
                    p::RigidParams,
                    material,
                    block::Int64)
    @info "Rigid material is applied. No internal forces are calculated. No deformation occurs only rigid body motion."
end

"""
    compute_model(nodes::AbstractVector{Int64}, p, material, time::Float64, dt::Float64)

Calculate the elastic bond force for each node.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p`: Model parameters; `material::BlockMaterial`: typed block material (`material.base`, `material.moduli`, `material.symmetry`).
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64},
                       p::RigidParams,
                       material,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
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

end
