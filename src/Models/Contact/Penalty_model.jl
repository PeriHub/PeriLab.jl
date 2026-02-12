# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Penalty_Model

using .....Data_Manager
using .....Helpers: get_shared_horizon, dot, norm

export contact_model_name
export init_contact_model
export compute_contact_model

function contact_model_name()
    return "Penalty Contact"
end

function init_contact_model(params)
    if !haskey(params, "Contact Stiffness")
        @warn "No ''Contact Stiffness'' has been defined. It is set to 1e8."
        params["Contact Stiffness"] = 1e8
    end
    @info "Contact Stiffness $(params["Contact Stiffness"])"
    if !haskey(params, "Friction Coefficient")
        params["Friction Coefficient"] = 0
    else
        if params["Friction Coefficient"] < 0
            @error "The Friction Coefficient must be greater or equal zero."
            return nothing
        end
    end

    if !haskey(params, "Symmetry")
        params["Symmetry"] = "3D"
    end
end
"""
    compute_contact_model(cg, params, compute_master_force_density,
                               compute_slave_force_density)
    Computes a Penalty model taken from [Peridigm](https://github.com/peridigm/peridigm/blob/master/src/contact/Peridigm_ShortRangeForceContactModel.cpp)

"""
function compute_contact_model(cg, params, compute_master_force_density,
                               compute_slave_force_density)
    contact_dict = Data_Manager.get_contact_dict(cg)
    contact_stiffness = params["Contact Stiffness"]
    contact_radius = params["Contact Radius"]
    normal_force = zeros(Data_Manager.get_dof())
    #Data_Manager.get_symmetry() # TODO store in materials the information
    for (master_id, contact) in pairs(contact_dict)
        for id in 1:contact["nSlaves"]
            slave_id = contact["Slaves"][id]
            horizon = get_shared_horizon(slave_id) # needed to get the correct contact horizon
            # TODO symmetry needed
            if params["Symmetry"] == "plane stress"
                stiffness = 9 / (pi * horizon^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
            elseif params["Symmetry"] == "plane strain"
                stiffness = 48 / (5 * pi * horizon^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
            else
                stiffness = 9 / (pi * horizon^5)  # -> from Peridigm
                #stiffness = 12 / (pi * horizon^4)  # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
            end

            @views distance = contact["Distances"][id]
            @views normal = contact["Normals"][id, :]
            temp = contact_stiffness * stiffness * (contact_radius - distance)
            normal_force = temp .* normal
            friction_id,
            friction_slave_id = compute_friction(id, slave_id,
                                                 params["Friction Coefficient"],
                                                 normal_force, normal)

            compute_master_force_density(master_id, slave_id,
                                         normal_force .+ friction_id)
            compute_slave_force_density(slave_id, master_id,
                                        normal_force .+ friction_slave_id)
        end
    end
end
"""

code was taken from Peridigm

"""
function compute_friction(id, slave_id, friction_coefficient, normal_force,
                          normal)
    dof = Data_Manager.get_dof()
    friction_id = zeros(dof)
    friction_slave_id = zeros(dof)
    if friction_coefficient == 0
        return friction_id, friction_slave_id
    end

    mapping = Data_Manager.get_exchange_id_to_local_id()

    if isnothing(get(mapping, id, nothing)) || isnothing(get(mapping, slave_id, nothing))
        return friction_id, friction_slave_id
    end

    velocity = Data_Manager.get_field("Velocity", "NP1")

    velo_id = zeros(dof)
    velo_slave_id = zeros(dof)
    vel_cm = zeros(dof)
    norm_vel_id = Float64(0)
    norm_velo_slave_id = Float64(0)

    @views current_dot_normal = dot(velocity[mapping[id], :], normal)
    @views current_dot_neighbor = dot(velocity[mapping[slave_id], :], normal)
    @views velo_id = velocity[mapping[id], :] .- current_dot_normal .* normal
    @views velo_slave_id = velocity[mapping[slave_id], :] .- current_dot_neighbor .* normal
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
