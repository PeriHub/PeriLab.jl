# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Zero_Energy_Control
using TimerOutputs: @timeit

using .....Data_Manager
using .....Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "control_name")
for mod in module_list
    include(mod["File"])
end

function init_model(nodes::AbstractVector{Int64}, material_parameter::Dict, block::Int64)
    if haskey(material_parameter, "Zero Energy Control")
        zero_energy_model = material_parameter["Zero Energy Control"]
        @info "Init zero energy control model ''$zero_energy_model'' at block $block."
        Data_Manager.set_analysis_model("Zero Energy Control Model", block,
                                        zero_energy_model)

        mod = create_module_specifics(zero_energy_model,
                                      module_list,
                                      @__MODULE__,
                                      "control_name")
        Data_Manager.set_model_module(zero_energy_model, mod)
        mod.init_model(nodes, material_parameter)
    else
        Data_Manager.set_analysis_model("Zero Energy Control Model", block, "")
        @warn "No zero energy control activated for corresponcence in block $block. This might cause errors."
    end
end

function compute_control(nodes::AbstractVector{Int64},
                         material_parameter::Dict{String,Any},
                         block::Int64,
                         time::Float64,
                         dt::Float64)
    for zero_energy_model in Data_Manager.get_analysis_model("Zero Energy Control Model",
                                        block)
        if zero_energy_model == ""
            continue
        end
        mod = Data_Manager.get_model_module(zero_energy_model)

        mod.compute_control(nodes,
                            material_parameter,
                            time,
                            dt)
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
