# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Zero_Energy_Control
using TimerOutputs: @timeit

using .....Data_Manager
using .....ModuleLoader: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "control_name")
for mod in module_list
    include(mod["File"])
end

"""
    init_model(nodes, material, block)

Initializes the zero energy control of a block (`material.base.zero_energy_control`).

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `block::Int64`: The block.
"""
function init_model(nodes::AbstractVector{Int64}, material, block::Int64)
    zero_energy_model = material.base.zero_energy_control
    if zero_energy_model !== nothing
        @debug "Init zero energy control model ''$zero_energy_model'' at block $block."
        Data_Manager.set_analysis_model("Zero Energy Control Model", block,
                                        zero_energy_model)
        mod = create_module_specifics(zero_energy_model,
                                      module_list,
                                      @__MODULE__,
                                      "control_name")
        Data_Manager.set_model_module(zero_energy_model, mod)
        mod.init_model(nodes, material)
    else
        Data_Manager.set_analysis_model("Zero Energy Control Model", block, "")
        @warn "No zero energy control activated for corresponcence in block $block. This might cause errors."
    end
end

"""
    compute_zero_energy_control(nodes, material, block, time, dt)

Applies the zero energy control of a block.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `block::Int64`: The block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_zero_energy_control(nodes::AbstractVector{Int64}, material, block::Int64,
                                     time::Float64, dt::Float64)
    for zero_energy_model in Data_Manager.get_analysis_model("Zero Energy Control Model",
                                                             block)
        zero_energy_model == "" && continue
        mod = Data_Manager.get_model_module(zero_energy_model)
        mod.compute_control(nodes, material, time, dt)
    end
end

```
create_zero_energy_mode_stiffness! interface for matrix based models

```
function create_zero_energy_mode_stiffness!(nodes::AbstractVector{Int64},
                                            dof::Int64,
                                            CVoigt::AbstractArray{Float64,3},
                                            Kinv::Array{Float64,3},
                                            zStiff::Array{Float64,3})
    zero_energy_model = Data_Manager.get_analysis_model("Zero Energy Control Model", 1)
    mod = Data_Manager.get_model_module(zero_energy_model[1])
    return mod.create_zero_energy_mode_stiffness!(nodes, dof, CVoigt, Kinv, zStiff)
end
end
