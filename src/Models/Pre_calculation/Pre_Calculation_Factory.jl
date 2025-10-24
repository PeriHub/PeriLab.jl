# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_Calculation

using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "pre_calculation_name")
for mod in module_list
    include(mod["File"])
end

using DataStructures

export init_fields
export init_model
export fields_for_local_synchronization
export compute_model
export check_dependencies

"""
    init_fields(datamanager::Module)

Initializes the fields.

# Arguments
- `datamanager::Data_Manager`: Datamanager
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function init_fields(datamanager::Module)
    dof = datamanager.get_dof()
    deformed_coorN,
    deformed_coorNP1 = datamanager.create_node_field("Deformed Coordinates",
                                                     Float64, dof)
    deformed_coorN = copy(datamanager.get_field("Coordinates"))
    deformed_coorNP1 = copy(datamanager.get_field("Coordinates"))
    datamanager.create_bond_field("Deformed Bond Geometry", Float64, dof)
    datamanager.create_bond_field("Deformed Bond Length", Float64, 1)
    datamanager.create_node_field("Displacements", Float64, dof)
    return datamanager
end

"""
    init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the model.

# Arguments
- `datamanager::Data_Manager`: Datamanager
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function init_model(datamanager::Module, nodes::AbstractVector{Int64},
                    block::Int64)
    dof = datamanager.get_dof()
    ## das muss hier rein. Das ist keine Komfortfunktion, sondern setzt Abhängigkeiten

    for (active_model_name, active_model) in pairs(datamanager.get_properties(block,
                                         "Pre Calculation Model"))
        if active_model
            mod = create_module_specifics(active_model_name,
                                          module_list,
                                          @__MODULE__,
                                          "pre_calculation_name")

            datamanager.set_model_module(active_model_name, mod)
            # TODO right now no additional information is needed
            # the parameter Dict is a kind of placeholder
            datamanager = mod.init_model(datamanager,
                                         nodes,
                                         Dict(datamanager.get_properties(block,
                                                                         "Pre Calculation Model")),
                                         block)
        end
    end

    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64)

Computes the pre calculation models

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       model_param::Union{Dict,OrderedDict},
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    for (pre_calculation_model, active) in pairs(model_param)
        if !active
            continue
        end
        mod = datamanager.get_model_module(pre_calculation_model)

        @timeit "compute $pre_calculation_model" datamanager=mod.compute(datamanager,
                                                                         nodes,
                                                                         model_param,
                                                                         block)
    end

    return datamanager
end

"""
    fields_for_local_synchronization(datamanager, model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `datamanager::Module`: datamanager.
- `model::String`: Model class.
- `block::Int64`: block ID
# Returns
- `datamanager::Module`: Datamanager.
"""

function fields_for_local_synchronization(datamanager, model, block)
    model_param = datamanager.get_properties(block, "Pre Calculation Model")

    for (pre_calculation_model, active) in pairs(model_param)
        if !active
            continue
        end
        mod = datamanager.get_model_module(pre_calculation_model)
        mod.fields_for_local_synchronization(datamanager, model)
    end
    return datamanager
end

"""
    check_dependencies(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}

Check if materials are used which needs a form of pre calculation. If so, the option will be set.

# Arguments
- `datamanager::Module`: Datamanager.
- `block_nodes::Dict{Int64,Vector{Int64}}`: block nodes.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function check_dependencies(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}})
    for block_id in eachindex(block_nodes)
        if !datamanager.check_property(block_id, "Material Model")
            continue
        end
        model_param = datamanager.get_properties(block_id, "Material Model")
        params_dict = datamanager.get_properties(block_id, "Pre Calculation Model")
        datamanager.set_properties(block_id,
                                   "Pre Calculation Model",
                                   merge(params_dict,
                                         Dict("Deformed Bond Geometry" => true)))
        if occursin("Correspondence", model_param["Material Model"])
            if haskey(model_param, "Bond Associated") && model_param["Bond Associated"]
                params_dict = datamanager.get_properties(block_id, "Pre Calculation Model")
                datamanager.set_properties(block_id,
                                           "Pre Calculation Model",
                                           merge(params_dict,
                                                 Dict("Bond Associated Correspondence" => true)))
                continue
            end
            params_dict = datamanager.get_properties(block_id, "Pre Calculation Model")
            datamanager.set_properties(block_id,
                                       "Pre Calculation Model",
                                       merge(params_dict,
                                             Dict("Shape Tensor" => true,
                                                  "Deformation Gradient" => true)))
        end
        # Check dependencies inside the pre calculation
        for (active_model_name, active_model) in pairs(datamanager.get_properties(block_id,
                                             "Pre Calculation Model"))
            if !active_model
                continue
            end
            if active_model_name == "Deformation Gradient"
                params_dict = datamanager.get_properties(block_id, "Pre Calculation Model")
                datamanager.set_properties(block_id,
                                           "Pre Calculation Model",
                                           merge(params_dict,
                                                 Dict("Deformed Bond Geometry" => true,
                                                      "Deformation Gradient" => true,
                                                      "Shape Tensor" => true)))
                continue
            end
            if active_model_name == "Shape Tensor"
                params_dict = datamanager.get_properties(block_id, "Pre Calculation Model")
                datamanager.set_properties(block_id,
                                           "Pre Calculation Model",
                                           merge(params_dict,
                                                 Dict("Deformed Bond Geometry" => true,
                                                      "Shape Tensor" => true)))
                continue
            end
            if active_model_name == "Bond Associated Correspondence"
                # Makes sure, that the neighbor information exists
                for local_block_id in eachindex(block_nodes)
                    params_dict = datamanager.get_properties(local_block_id,
                                                             "Pre Calculation Model")
                    datamanager.set_properties(block_id,
                                               "Pre Calculation Model",
                                               merge(params_dict,
                                                     Dict("Deformed Bond Geometry" => true,
                                                          "Bond Associated Correspondence" => true)))
                end
                continue
            end
        end
        # sort all elements in the pre defined order
        order_dict = datamanager.get_properties(block_id, "Pre Calculation Model")
        order_vector = datamanager.get_pre_calculation_order()
        datamanager.set_properties(block_id,
                                   "Pre Calculation Model",
                                   OrderedDict(k => order_dict[k]
                                               for k in order_vector
                                               if haskey(order_dict, k)))
    end
    return datamanager
end

end
