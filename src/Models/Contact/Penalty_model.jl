# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Penalty_model

include("../../Support/Helpers.jl")
import .Helpers: get_shared_horizon, dot, norm

export contact_model_name
export init_contact_model
export compute_contact_model

function contact_model_name()
    return "Penalty Contact"
end

function init_contact_model(datamanager, params)
    if !haskey(params, "Contact Stiffness")
        @warn "No ''Contact Stiffness'' has been defined. It is set to 1e8."
        params["Contact Stiffness"] = 1e8
    end
    if !haskey(params, "Friction Coefficient")
        params["Friction Coefficient"] = 0
    else
        if params["Friction Coefficient"] < 0
            @error "The Friction Coefficient must be greater or equal zero."
            return nothing
        end
    end
    return datamanager
end

function compute_contact_model(datamanager, cg, params, compute_master_force_density,
                               compute_slave_force_density)
    contact_dict = datamanager.get_contact_dict(cg)
    contact_stiffness = params["Contact Stiffness"]
    contact_radius = params["Contact Radius"]
    normal_force = zeros(datamanager.get_dof())
    #datamanager.get_symmetry() # TODO store in materials the information
    for (master_id, contact) in pairs(contact_dict)
        for id in 1:contact["nSlaves"]
            slave_id = contact["Slaves"][id]
            horizon = get_shared_horizon(datamanager, slave_id) # needed to get the correct contact horizon
            # TODO symmetry needed
            stiffness = 48.0 / 5.0 * contact_stiffness / (pi * horizon^3)

            @views distance = contact["Distances"][id]
            @views normal = contact["Normals"][id, :]
            temp = (contact_radius - distance) / horizon
            normal_force = stiffness * temp .* normal
            friction_id,
            friction_slave_id = compute_friction(datamanager, id, slave_id,
                                                 params["Friction Coefficient"],
                                                 normal_force, normal)

            compute_master_force_density(datamanager, master_id, slave_id,
                                         normal_force .+ friction_id)
            compute_slave_force_density(datamanager, slave_id, master_id,
                                        normal_force .+ friction_slave_id)
        end
    end
    return datamanager
end

function compute_friction(datamanager, id, slave_id, friction_coefficient, normal_force,
                          normal)
    dof = datamanager.get_dof()
    friction_id = zeros(dof)
    friction_slave_id = zeros(dof)
    if friction_coefficient == 0
        return friction_id, friction_slave_id
    end

    mapping = datamanager.get_exchange_id_to_local_id()

    if !(id in keys(mapping)) || !(slave_id in keys(mapping))
        return
    end

    velocity = datamanager.get_field("Velocity", "NP1")

    velo_id = zeros(dof)
    velo_slave_id = zeros(dof)
    vel_cm = zeros(dof)
    norm_vel_id = Float64(0)
    norm_velo_slave_id = Float64(0)

    @views current_dot_normal = dot(velocity[mapping[id], :], normal)
    @views current_dot_neighbor = dot(velocity[mapping[slave_id], :], normal)
    @views velo_id = velocity[mapping[id], :] .- current_dot_normal .* normal
    @views velo_slave_id = velocity[mapping[slave_id], :] .- current_dot_normal .* normal
    @views vel_cm = 0.5 .* (velo_slave_id + velo_id)
    @views velo_id .-= vel_cm
    @views velo_slave_id .-= vel_cm
    norm_vel_id = norm(velo_id)
    norm_velo_slave_id = norm(velo_slave_id)

    if norm_vel_id != 0
        @views friction_id = -friction_coefficient * norm(normal_force) .* velo_id ./
                             norm_vel_id
    end
    if norm_velo_slave_id != 0
        @views friction_slave_id = -friction_coefficient * norm(normal_force) .*
                                   velo_slave_id ./ norm_velo_slave_id
    end
    return friction_id, friction_slave_id
end
end
