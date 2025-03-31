# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact
include("Contact_search.jl")
using .Contact_search: init_contact_search
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "contact_model_name")
Set_modules.include_files(module_list)
export init_model
export compute_model

"""
    init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the contact model.

# Arguments
- `datamanager::Data_manager`: Datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module, params, nodes::Union{SubArray,Vector{Int64}})
    surface_nodes = fill(Vector{Int64}([]), length(params))
    for (cm, contact_params) in enumerate(params)
        init_contact_search(datamanager, contact_params, cm, surface_nodes)
    end

    #datamanager.create_constant_field("Contact IDs", Int64, ??)
    # allnodes
    # cm -> check if inside or not
    #  reduce allnodes

    # surface_nodes
    # global ids; all nodes reduced global_id => local_id <=> contact_global
    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, block::Int64, time::Float64, dt::Float64, to::TimerOutput)

Compute the forces of the contact model.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `model_param::Dict`: The contact parameter.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_model(datamanager::Module,
                       nodes::Union{SubArray,Vector{Int64}},
                       model_param::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    return datamanager
end
end
