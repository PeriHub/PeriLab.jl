# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal

using ....Data_Manager
using ....PeriLabExceptions: @abort
using TimerOutputs: @timeit
using ....ModuleLoader: find_registered_modules, licensed_modules
using ....Helpers: each_call

using .....ParameterSpec: @params, register_base!, WithBase, model_parts, model_module,
                          model_function, model_label

"""
    ThermalBaseParams

Keys every thermal model may use: `Thermal Conductivity` (read by the critical
time step for any thermal model, and by Thermal Flow).
"""
@params "Thermal model of the blocks that use this entry, with the keys shared by all thermal models" struct ThermalBaseParams
    "Thermal conductivity of the material; required by Thermal Flow and used for the stable time step"
    thermal_conductivity::Union{Nothing,Float64} = opt("Thermal Conductivity";
                                                       default = nothing, min = 0)
end

# registration runs at load time, never during precompilation
__init__() = register_base!(:thermal, ThermalBaseParams)
for file in find_registered_modules(@__DIR__, "register_thermal")
    include(file)
end

global licensed = Any[]

"""
    load_licensed_models()

Includes the licensed `Interlayer` helper and thermal modules (if a license is configured)
so that they register their parameter structs before the input deck is read. Runs once.
"""
function load_licensed_models()
    isempty(licensed) || return licensed
    global licensed = vcat(licensed_modules(@__MODULE__, "Interlayer"),
                           licensed_modules(@__MODULE__, "Thermal"))
    return licensed
end

"""
    load_licensed_models(raw_deck)

Loads the licensed thermal modules only if the raw input deck (the YAML dict) has
`Models: Thermal Models`, so that other runs never contact the license server.
"""
function load_licensed_models(raw_deck::AbstractDict)
    models = get(get(raw_deck, "PeriLab", Dict()), "Models", nothing)
    (models isa AbstractDict && haskey(models, "Thermal Models")) || return Any[]
    return load_licensed_models()
end

# not a deck (empty file, a list): the input validation reports it
load_licensed_models(raw_deck) = Any[]

"""
    post_process_licensed(input, exodus_files::Vector{String})

Calls `post_process(input, exodus_files)` of every loaded licensed thermal module that
defines it, after the simulation (rank 0, final Exodus files).
"""
function post_process_licensed(input, exodus_files::Vector{String})
    # invokelatest: the licensed modules are defined at runtime (newer world)
    Base.invokelatest() do
        for m in licensed
            mod = getfield(@__MODULE__, Symbol(m["Module Name"]))
            isdefined(mod, :post_process) && mod.post_process(input, exodus_files)
        end
    end
end

export init_model
export compute_model
export init_fields
export fields_for_local_synchronization

"""
    init_fields()

Initialize thermal model fields
"""
function init_fields()
    Data_Manager.create_node_scalar_field("Temperature", Float64)
    Data_Manager.create_constant_node_scalar_field("Delta Temperature", Float64)
    Data_Manager.create_node_scalar_field("Heat Flow", Float64)
    # ghost nodes are not integrated on this rank; without it they stay at 0 K and act as heat sinks
    Data_Manager.set_synch("Temperature", false, true)
    Data_Manager.create_constant_node_scalar_field("Specific Volume", Float64)
    # if it is already initialized via mesh file no new field is created here
    Data_Manager.create_constant_node_scalar_field("Surface_Nodes", Bool;
                                                   default_value = true)
end

"""
    compute_model(nodes::AbstractVector{Int64}, thermal::WithBase, block::Int64, time::Float64, dt::Float64)

Computes every part of the block's thermal model, in deck order.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `thermal::WithBase`: The typed thermal model of the block (`thermal.base`: Thermal Conductivity)
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64}, thermal::WithBase, block::Int64,
                       time::Float64, dt::Float64)
    # the parts of a composite have different types: a loop would box them
    each_call(compute_part!, model_parts(thermal.model), nodes, thermal, block, time, dt)
    return nothing
end

function compute_part!(part, nodes::AbstractVector{Int64}, thermal::WithBase, block::Int64,
                       time::Float64, dt::Float64)
    @timeit model_label(part) model_function(part, Val(:compute_model))(nodes, part,
                                                                        thermal, block,
                                                                        time, dt)
    return nothing
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initializes every part of the block's thermal model (`Data_Manager.get_block_models(block).thermal`).

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    thermal = Data_Manager.get_block_models(block).thermal
    for part in model_parts(thermal.model)
        # invokelatest: licensed models are defined at runtime (compute_model runs inside
        # the solver's invokelatest already)
        Base.@invokelatest model_module(part).init_model(nodes, part, thermal, block)
    end
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block ID
"""
function fields_for_local_synchronization(model, block)
    thermal = Data_Manager.get_block_models(block).thermal
    for part in model_parts(thermal.model)
        Base.@invokelatest model_module(part).fields_for_local_synchronization(model)
    end
end

end
