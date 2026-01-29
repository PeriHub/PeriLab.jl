# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Model_Factory

using TimerOutputs: @timeit
using ...Data_Manager
using ...Helpers:
                  check_inf_or_nan, find_active_nodes, get_active_update_nodes, invert,
                  determinant, matrix_style, eigvals
include("./Pre_calculation/Pre_Calculation_Factory.jl")
include("./Surface_correction/Surface_correction.jl")
include("./Contact/Contact_Factory.jl")
include("./Additive/Additive_Factory.jl")
include("./Degradation/Degradation_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
using ...Parameter_Handling: get_model_parameter, get_heat_capacity
using .Additive
using .Degradation
using .Damage
using .Material
using .Pre_Calculation
using .Surface_Correction: init_surface_correction, compute_surface_correction
using .Thermal
using .Contact
# in future FEM will be outside of the Model_Factory
using ..FEM
export compute_models
export init_models
export read_properties

"""
	init_models(params::Dict, block_nodes::Dict{Int64,Vector{Int64}}, solver_options::Dict)

Initialize models

# Arguments
- `params::Dict`: Parameters.
- `block_nodes::Dict{Int64,Vector{Int64}}`: block nodes.
- `solver_options::Dict`: Solver options.
"""
function init_models(params::Dict,
                     block_nodes::Dict{Int64,Vector{Int64}},
                     solver_options::Dict,
                     synchronise_field)
    if "Pre_Calculation" in solver_options["Models"]
        @info "Check pre calculation models are initialized for material models"
        Pre_Calculation.check_dependencies(block_nodes)
        if haskey(params["Models"], "Material Models")
            for mat in keys(params["Models"]["Material Models"])
                if haskey(params["Models"]["Material Models"][mat], "Accuracy Order")
                    Data_Manager.set_accuracy_order(params["Models"]["Material Models"][mat]["Accuracy Order"])
                end
            end
        end
    end
    for model_name in solver_options["Models"]
        add_model(model_name)
    end
    for model_name in solver_options["All Models"]
        add_model(model_name, true)
    end

    if "Additive" in solver_options["Models"] || "Thermal" in solver_options["Models"]
        heat_capacity = Data_Manager.create_constant_node_scalar_field("Specific Heat Capacity",
                                                                       Float64)
        heat_capacity = set_heat_capacity(params, block_nodes, heat_capacity) # includes the neighbors
    end

    if isnothing(Data_Manager.get_step()) || Data_Manager.get_step() == 1
        for (active_model_name, active_model) in pairs(Data_Manager.get_active_models(true))
            @debug "Init $active_model_name fields"
            @timeit "$active_model_name model fields" active_model.init_fields()
        end
    end
    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())
        @info "Init $active_model_name"

        model_used = false
        for block in eachindex(block_nodes)
            if Data_Manager.check_property(block, active_model_name)
                @timeit "init $active_model_name models" active_model.init_model(block_nodes[block],
                                                                                 block)
                @timeit "init fields_for_local_synchronization $active_model_name models" active_model.fields_for_local_synchronization(active_model_name,
                                                                                                                                        block)
                if active_model_name == "Damage Model" &&
                   haskey(Data_Manager.get_properties(block, active_model_name),
                          "Local Damping")
                    Material.init_local_damping(block_nodes[block],
                                                Data_Manager.get_properties(block,
                                                                            "Material Model"),
                                                Data_Manager.get_properties(block,
                                                                            "Damage Model"))
                end
                model_used = true
                # put it in Data_Manager
            end
        end
        if !model_used
            @warn "$active_model_name is defined and activated, but not referenced in any block definition"
        end
    end

    init_surface_correction(params, local_synch,
                            synchronise_field)

    if solver_options["Calculation"]["Calculate Cauchy"] |
       solver_options["Calculation"]["Calculate von Mises stress"]
        Data_Manager.create_node_tensor_field("Cauchy Stress",
                                              Float64, Data_Manager.get_dof())
    end
    if solver_options["Calculation"]["Calculate Strain"]
        Data_Manager.create_node_tensor_field("Strain", Float64, Data_Manager.get_dof())
    end
    if solver_options["Calculation"]["Calculate von Mises stress"]
        Data_Manager.create_node_scalar_field("von Mises Stress", Float64)
    end

    check_contact(params)
    @info "Finalize Init Models"
end

function check_contact(params::Dict)
    if haskey(params, "Contact")
        return Contact.init_contact_model(params["Contact"])
    end
end
function check_contact(params::Dict, time::Float64, dt::Float64)
    if length(params) != 0
        return Contact.compute_contact_model(params, time,
                                             dt)
    end
end

"""
	compute_models(block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Vector{String}, synchronise_field)

Computes the material point models

# Arguments
- `block_nodes::Dict{Int64,Vector{Int64}}`: The block nodes
- `dt::Float64`: The time step
- `time::Float64`: The current time of the solver
- `options::Vector{String}`: The options
- `synchronise_field`: The synchronise field
"""
function compute_models(block_nodes::Dict{Int64,Vector{Int64}},
                        dt::Float64,
                        time::Float64,
                        options::Vector{String},
                        synchronise_field)
    fem_option = Data_Manager.fem_active()
    if fem_option
        fe_nodes = Data_Manager.get_field("FE Nodes")
    end

    active_list = Data_Manager.get_field("Active")
    # TODO check if pre calculation should run block wise. For mixed model applications it makes sense.
    # TODO add for pre calculation a whole model option, to get the neighbors as well, e.g. for bond associated
    # TODO check for loop order?
    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())
        #synchronise_field(Data_Manager.local_synch_fiels(active_model_name))
        if active_model_name == "Damage Model"
            continue
        end

        local_synch(active_model_name, "upload_to_cores", synchronise_field)
        # maybe not needed?
        local_synch(active_model_name,
                    "download_from_cores",
                    synchronise_field)

        for (block, nodes) in pairs(block_nodes)
            # "delete" the view of active nodes
            active_nodes = Data_Manager.get_field("Active Nodes")

            active_nodes = find_active_nodes(active_list,
                                             active_nodes,
                                             nodes,
                                             active_model_name != "Additive Model")

            if fem_option # TODO might lead to problems in 3D
                # find all non-FEM nodes in active nodes
                active_nodes = Data_Manager.get_field("Active Nodes")
                active_nodes = find_active_nodes(fe_nodes,
                                                 active_nodes,
                                                 1:Data_Manager.get_nnodes(),
                                                 false)
                if active_nodes == []
                    continue
                end
            end
            if Data_Manager.check_property(block, active_model_name)
                # synch
                @timeit "compute $active_model_name" active_model.compute_model(active_nodes,
                                                                                Data_Manager.get_properties(block,
                                                                                                            active_model_name),
                                                                                block,
                                                                                time,
                                                                                dt)
            end
        end
    end

    # Why not update_list.=false? -> avoid neighbors
    update_list = Data_Manager.get_field("Update")
    for (block, nodes) in pairs(block_nodes)
        update_list[nodes] .= false
    end

    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())
        if active_model_name == "Additive Model"
            continue
        end
        local_synch(active_model_name, "upload_to_cores", synchronise_field)
        # maybe not needed?
        local_synch(active_model_name,
                    "download_from_cores",
                    synchronise_field)
        for (block, nodes) in pairs(block_nodes)
            active_nodes = Data_Manager.get_field("Active Nodes")
            update_nodes = Data_Manager.get_field("Update Nodes")
            active_nodes = find_active_nodes(active_list, active_nodes, nodes)
            if fem_option
                # FEM active means FEM nodes
                active_nodes = Data_Manager.get_field("Active Nodes")
                active_nodes = find_active_nodes(fe_nodes,
                                                 active_nodes,
                                                 1:Data_Manager.get_nnodes(),
                                                 false)
                if active_nodes == []
                    continue
                end
            end
            update_nodes = get_update_nodes(active_list,
                                            update_list,
                                            nodes,
                                            update_nodes,
                                            active_nodes,
                                            active_model_name)

            # active or all, or does it not matter?

            if Data_Manager.check_property(block, active_model_name)
                # TODO synch
                @timeit "compute $active_model_name" active_model.compute_model(update_nodes,
                                                                                Data_Manager.get_properties(block,
                                                                                                            active_model_name),
                                                                                block,
                                                                                time,
                                                                                dt)
            end
        end
    end
    # must be here to avoid double distributions
    # distributes ones over all nodes
    if fem_option
        @timeit "FEM" begin
            nelements = Data_Manager.get_num_elements()

            @timeit "eval" FEM.eval_FEM(Vector{Int64}(1:nelements),
                                        Data_Manager.get_properties(1, "FEM"),
                                        time,
                                        dt)
            active_nodes = Data_Manager.get_field("Active Nodes")

            FEM.force_densities(find_active_nodes(fe_nodes, active_nodes,
                                                  1:Data_Manager.get_nnodes(), true))
        end
    end

    if "Material" in options
        if "Damage" in options
            for (block, nodes) in pairs(block_nodes)
                if haskey(Data_Manager.get_properties(block, "Damage Model"),
                          "Local Damping")
                    active_nodes = Data_Manager.get_field("Active Nodes")
                    if fem_option
                        active_nodes = find_active_nodes(active_list,
                                                         active_nodes,
                                                         find_active_nodes(fe_nodes,
                                                                           active_nodes,
                                                                           nodes))
                    else
                        find_active_nodes(active_list, active_nodes, nodes)
                    end
                    @timeit "local_damping_due_to_damage" Material.compute_local_damping(active_nodes,
                                                                                         Data_Manager.get_properties(block,
                                                                                                                     "Damage Model")["Local Damping"],
                                                                                         dt)
                end
            end
        end
        active_nodes = Data_Manager.get_field("Active Nodes")
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:Data_Manager.get_nnodes())

        compute_surface_correction(active_nodes,
                                   local_synch,
                                   synchronise_field)

        @timeit "distribute_force_densities" Material.distribute_force_densities(active_nodes)
    end

    if fem_option
        @timeit "coupling" FEM.Coupling.compute_coupling(Data_Manager.get_properties(1,
                                                                                     "FEM"))
    end
    check_contact(Data_Manager.get_contact_properties(), time, dt)
    #=
    Used for shape tensor or other fixed calculations, to avoid an update if its not needed.
    The damage update is done in the second loop.
    =#
    update_list = Data_Manager.get_field("Update")
    for (block, nodes) in pairs(block_nodes)
        update_list[nodes] .= false
    end
end

"""
	compute_stiff_matrix_compatible_models(block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Vector{String}, synchronise_field)

Computes the models models that are compatible with the stiffness matrix calculation.

# Arguments
- `block_nodes::Dict{Int64,Vector{Int64}}`: The block nodes
- `dt::Float64`: The time step
- `time::Float64`: The current time of the solver
- `options::Vector{String}`: The options
- `synchronise_field`: The synchronise field
"""
function compute_stiff_matrix_compatible_models(block_nodes::Dict{Int64,Vector{Int64}},
                                                dt::Float64,
                                                time::Float64,
                                                options::Vector{String},
                                                synchronise_field)
    active_list = Data_Manager.get_field("Active")

    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())

        #local_synch(Data_Manager, active_model_name, "upload_to_cores", synchronise_field)
        # maybe not needed?
        #local_synch(Data_Manager,
        #	active_model_name,
        #	"download_from_cores",
        #	synchronise_field)
        #TODO: Change Thermal for Thermal Expansion only!
        if active_model_name == "Material Model" && !("Thermal" in options)
            # we need here an activation trigger for mixed models in future
            continue
        end

        for (block, nodes) in pairs(block_nodes)
            # "delete" the view of active nodes
            active_nodes = Data_Manager.get_field("Active Nodes")

            active_nodes = find_active_nodes(active_list,
                                             active_nodes,
                                             nodes,
                                             active_model_name != "Additive Model")

            if Data_Manager.check_property(block, active_model_name)
                # synch
                @timeit "compute $active_model_name" active_model.compute_model(active_nodes,
                                                                                Data_Manager.get_properties(block,
                                                                                                            active_model_name),
                                                                                block,
                                                                                time,
                                                                                dt)
            end
        end
    end

    if ("Material" in options) && ("Thermal" in options)
        active_nodes = Data_Manager.get_field("Active Nodes")
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:Data_Manager.get_nnodes())
        @timeit "distribute_force_densities" Material.distribute_force_densities(active_nodes)
    end
end

function get_update_nodes(active_list,
                          update_list,
                          nodes,
                          update_nodes,
                          active_nodes,
                          active_model_name)
    if active_model_name == "Damage Model"
        return @view active_nodes[:]
    else
        return get_active_update_nodes(active_list, update_list, nodes, update_nodes)
    end
end

"""
	get_block_model_definition(params::Dict, block_id_list::Int64, prop_keys::Vector{String}, properties)

Get block model definition.

Special case for pre calculation. It is set to all blocks, if no block definition is defined, but pre calculation is.

# Arguments
- `params::Dict`: Parameters.
- `block_id_list::Vector{Int64}`: List of block id's.
- `prop_keys::Vector{String}`: Property keys.
- `properties`: Properties function.
# Returns
- `properties`: Properties function.
"""
function get_block_model_definition(params::Dict,
                                    block_name_list::Vector{String},
                                    block_id_list::Vector{Int64},
                                    prop_keys::Vector{String},
                                    properties,
                                    directory::String = "",
                                    material_model::Bool = true)
    # properties function from Data_Manager

    if haskey(params["Models"], "Pre Calculation Global")
        for block_id in block_id_list
            properties(block_id,
                       "Pre Calculation Model",
                       params["Models"]["Pre Calculation Global"])
        end
    end

    for (block_id, block_name) in zip(block_id_list, block_name_list)
        if !haskey(params["Blocks"], block_name)
            continue
        end
        block = params["Blocks"][block_name]
        for model in prop_keys
            if model == "Material Model" && !material_model
                continue
            end
            if haskey(block, model)
                properties(block_id,
                           model,
                           get_model_parameter(params, model, block[model], directory))
            end
        end
    end
    return properties
end

"""
	read_properties(params::Dict, material_model::Bool)

Read properties of material.

# Arguments
- `params::Dict`: Parameters.
- `material_model::Bool`: Material model.
"""
function read_properties(params::Dict, material_model::Bool)
    Data_Manager.init_properties()
    block_name_list = Data_Manager.get_block_name_list()
    block_id_list = Data_Manager.get_block_id_list()
    prop_keys = Data_Manager.init_properties()
    directory = Data_Manager.get_directory()
    get_block_model_definition(params,
                               block_name_list,
                               block_id_list,
                               prop_keys,
                               Data_Manager.set_properties,
                               directory,
                               material_model)
    if material_model
        dof = Data_Manager.get_dof()
        for block in block_id_list
            Material.check_material_symmetry(dof,
                                             Data_Manager.get_properties(block,
                                                                         "Material Model"))
            Material.determine_isotropic_parameter(Data_Manager.get_properties(block,
                                                                               "Material Model"))
        end
    end
end

"""
	set_heat_capacity(params::Dict, block_nodes::Dict, heat_capacity::NodeScalarField{Float64})

Sets the heat capacity of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: The block nodes
- `heat_capacity::NodeScalarField{Float64}`: The heat capacity array
# Returns
- `heat_capacity::SubArray`: The heat capacity array
"""
function set_heat_capacity(params::Dict, block_nodes::Dict,
                           heat_capacity::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        heat_capacity[block_nodes[block]] .= get_heat_capacity(params, block)
    end
    return heat_capacity
end

"""
	add_model(model_name::String)

Includes the models in the Data_Manager and checks if the model definition is correct or not.

# Arguments
- `model_name::String`: The block nodes
"""
function add_model(model_name::String, all::Bool = false)
    try
        # to catch "Pre_Calculation"
        Data_Manager.add_active_model(replace(model_name, "_" => " ") * " Model",
                                      eval(Meta.parse(model_name)),
                                      all)
    catch
        @error "Model $model_name is not specified and cannot be included."
        return nothing
    end
end

function local_synch(model, direction, synchronise_field)
    synch_fields = Data_Manager.get_local_synch_fields(model)
    for synch_field in keys(synch_fields)
        synchronise_field(Data_Manager.get_comm(),
                          synch_fields,
                          Data_Manager.get_overlap_map(),
                          Data_Manager.get_field,
                          synch_field,
                          direction)
    end
end

"""
	compute_thermodynamic_critical_time_step(nodes::AbstractVector{Int64}, lambda::Float64, Cv::Float64)

Calculate the critical time step for a thermodynamic simulation based on  [OterkusS2014](@cite).

This function iterates over a collection of nodes and computes the critical time step for each node using provided input data and parameters.

# Arguments
- `nodes::AbstractVector{Int64}`: The collection of nodes to calculate the critical time step for.
- `lambda::Float64`: The material parameter used in the calculations.
- `Cv::Float64`: The heat capacity at constant volume used in the calculations.

# Returns
- `Float64`: The calculated critical time step for the thermodynamic simulation.

# Dependencies
This function depends on the following data fields from the `Data_Manager` module:
- `get_nlist()`: Returns the neighbor list.
- `get_field("Density")`: Returns the density field.
- `get_field("Bond Length")`: Returns the bond distance field.
- `get_field("Volume")`: Returns the volume field.
- `get_field("Number of Neighbors")`: Returns the number of neighbors field.
"""
function compute_thermodynamic_critical_time_step(nodes::AbstractVector{Int64},
                                                  lambda::Union{Float64,Int64})
    critical_time_step::Float64 = 1.0e50
    nlist = Data_Manager.get_nlist()
    density = Data_Manager.get_field("Density")
    undeformed_bond_length = Data_Manager.get_field("Bond Length")
    volume = Data_Manager.get_field("Volume")
    Cv = Data_Manager.get_field("Specific Heat Capacity")
    lambda = matrix_style(lambda)
    eigLam = maximum(eigvals(lambda))

    for iID in nodes
        denominator = get_cs_denominator(volume[nlist[iID]], undeformed_bond_length[iID])
        t = density[iID] * Cv[iID] / (eigLam * denominator)
        critical_time_step = test_timestep(t, critical_time_step)
    end
    return sqrt(critical_time_step)
end

"""
	get_cs_denominator(volume::AbstractVector{Float64}, undeformed_bond::AbstractVector{Float64})

Calculate the denominator for the critical time step calculation.

# Arguments
- `volume::AbstractVector{Float64}`: The volume field.
- `undeformed_bond::Union{SubArray,Vector{Float64},Vector{Int64}}`: The undeformed bond field.
# Returns
- `Float64`: The denominator for the critical time step calculation.
"""
function get_cs_denominator(volume::AbstractVector{Float64},
                            undeformed_bond::AbstractVector{Float64})::Float64
    return sum(volume ./ undeformed_bond)
end

"""
	compute_mechanical_critical_time_step(nodes::AbstractVector{Int64}, bulk_modulus::Float64)

Calculate the critical time step for a mechanical simulation using a bond-based approximation [LittlewoodDJ2013](@cite).

This function iterates over a collection of nodes and computes the critical time step for each node based on the given input data and parameters.

# Arguments
- `nodes::AbstractVector{Int64}`: The collection of nodes to calculate the critical time step for.
- `bulk_modulus::Float64`: The bulk modulus used in the calculations.

# Returns
- `Float64`: The calculated critical time step for the mechanical simulation.

# Dependencies
This function depends on the following data fields from the `Data_Manager` module:
- `get_nlist()`: Returns the neighbor list.
- `get_field("Density")`: Returns the density field.
- `get_field("Bond Length")`: Returns the bond distance field.
- `get_field("Volume")`: Returns the volume field.
- `get_field("Horizon")`: Returns the horizon field.
"""
function compute_mechanical_critical_time_step(nodes::AbstractVector{Int64},
                                               bulk_modulus::Union{Float64,Int64,SubArray,
                                                                   Vector{Float64}})
    critical_time_step::Float64 = 1.0e50
    nlist = Data_Manager.get_nlist()
    density = Data_Manager.get_field("Density")
    undeformed_bond_length = Data_Manager.get_field("Bond Length")
    volume = Data_Manager.get_field("Volume")
    horizon = Data_Manager.get_field("Horizon")

    for iID in nodes
        denominator = get_cs_denominator(volume[nlist[iID]], undeformed_bond_length[iID])
        # TODO Adapt to 2D applications
        springConstant = 18.0 * maximum(bulk_modulus) /
                         (pi * horizon[iID] * horizon[iID] * horizon[iID] * horizon[iID])

        t = density[iID] / (denominator * springConstant)
        critical_time_step = test_timestep(t, critical_time_step)
    end
    return sqrt(2 * critical_time_step)
end

"""
	test_timestep(t::Float64, critical_time_step::Float64)

Compare a time step `t` with a critical time step `critical_time_step` and update `critical_time_step` if `t` is smaller.

# Arguments
- `t::Float64`: The time step to compare with `critical_time_step`.
- `critical_time_step::Float64`: The current critical time step.

# Returns
- `critical_time_step::Float64`: The updated critical time step, which is either the original `critical_time_step` or `t`, whichever is smaller.
"""
function test_timestep(t::Float64, critical_time_step::Float64)
    if t < critical_time_step
        critical_time_step = t
    end
    return critical_time_step
end

"""
	compute_crititical_time_step(block_nodes::Dict{Int64,Vector{Int64}}, mechanical::Bool, thermo::Bool)

Calculate the critical time step for a simulation considering both mechanical and thermodynamic aspects.

This function computes the critical time step by considering mechanical and thermodynamic properties of different blocks. The resulting critical time step is based on the smallest critical time step found among the blocks.

# Arguments
- `block_nodes::Dict{Int64, Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `mechanical::Bool`: If `true`, mechanical properties are considered in the calculation.
- `thermo::Bool`: If `true`, thermodynamic properties are considered in the calculation.

# Returns
- `Float64`: The calculated critical time step based on the smallest critical time step found among the blocks.

# Dependencies
This function may depend on the following functions:
- `compute_thermodynamic_critical_time_step`: Used if `thermo` is `true` to calculate thermodynamic critical time steps.
- `compute_mechanical_critical_time_step`: Used if `mechanical` is `true` to calculate mechanical critical time steps.
- The availability of specific properties from the data manager module.

# Errors
- If required properties are not available in the data manager, it may raise an error message.
"""
function compute_crititical_time_step(block_nodes::Dict{Int64,Vector{Int64}},
                                      mechanical::Bool,
                                      thermal::Bool)
    critical_time_step::Float64 = 1.0e50
    for iblock in eachindex(block_nodes)
        if thermal
            lambda = Data_Manager.get_property(iblock, "Thermal Model",
                                               "Thermal Conductivity")
            # if Cv and lambda are not defined it is valid, because an analysis can take place, if material is still analysed
            if isnothing(lambda)
                if !mechanical
                    @warn "No time step can be calculated, because the heat conduction is not defined."
                    return critical_time_step
                end
            else
                t = compute_thermodynamic_critical_time_step(block_nodes[iblock],
                                                             lambda)
                critical_time_step = test_timestep(t, critical_time_step)
            end
        end
        if mechanical
            bulk_modulus = Data_Manager.get_property(iblock, "Material Model",
                                                     "Bulk Modulus")
            nu_xy = Data_Manager.get_property(iblock, "Material Model",
                                              "Poisson's Ratio XY")
            nu_yz = Data_Manager.get_property(iblock, "Material Model",
                                              "Poisson's Ratio YZ")
            nu_xz = Data_Manager.get_property(iblock, "Material Model",
                                              "Poisson's Ratio XZ")
            E_x = Data_Manager.get_property(iblock, "Material Model", "Young's Modulus X")
            E_y = Data_Manager.get_property(iblock, "Material Model", "Young's Modulus Y")
            E_z = Data_Manager.get_property(iblock, "Material Model", "Young's Modulus Z")
            g_xy = Data_Manager.get_property(iblock, "Material Model", "Shear Modulus XY")
            g_yz = Data_Manager.get_property(iblock, "Material Model", "Shear Modulus YZ")
            c_44 = Data_Manager.get_property(iblock, "Material Model", "C44")
            c_55 = Data_Manager.get_property(iblock, "Material Model", "C55")
            c_66 = Data_Manager.get_property(iblock, "Material Model", "C66")
            if !isnothing(bulk_modulus)
                bulk_modulus = bulk_modulus
            elseif !isnothing(nu_xy) && !isnothing(nu_yz) && !isnothing(nu_xz)
                s11 = 1 / E_x
                s22 = 1 / E_y
                s33 = 1 / E_z
                s12 = -nu_xy / E_x
                s23 = -nu_yz / E_z
                s13 = -nu_xz / E_z
                bulk_modulus = 1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s13))
            elseif !isnothing(c_44) && !isnothing(c_55) && !isnothing(c_66)
                bulk_modulus = maximum([c_44 / 2, c_55 / 2, c_66 / 2])
                #TODO: temporary solution!!!
            elseif !isnothing(g_xy)
                bulk_modulus = g_xy / 2
                #TODO: temporary solution!!!
            else
                @error "No time step for material is determined because of missing properties."
                return nothing
            end
            t = compute_mechanical_critical_time_step(block_nodes[iblock],
                                                      bulk_modulus)
            critical_time_step = test_timestep(t, critical_time_step)
        end
    end
    return critical_time_step
end

end
