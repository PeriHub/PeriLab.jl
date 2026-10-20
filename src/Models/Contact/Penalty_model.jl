# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Penalty_Model

using .....Data_Manager
using .....PeriLabExceptions: @abort
using .....Helpers: get_shared_horizon, dot, norm, with_dof
using .....ParameterSpec: @params, register_contact
using StaticArrays: SVector

export init_contact_model
export compute_contact_model

@params struct PenaltyContactParams
    contact_stiffness::Float64 = opt("Contact Stiffness"; default = 1e8, min = 0)
    friction_coefficient::Float64 = opt("Friction Coefficient"; default = 0.0, min = 0)
end

__init__() = register_contact("Penalty Contact", PenaltyContactParams)

"""
    init_contact_model(p, contact)

Initializes the penalty contact model. `contact` holds the shared contact keys
in `contact.base`.
"""
function init_contact_model(p::PenaltyContactParams, contact)
    @info "Contact Stiffness $(p.contact_stiffness)"
    return nothing
end
"""
    compute_contact_model(cg, p, contact, compute_master_force_density,
                               compute_slave_force_density)
    Computes a Penalty model taken from [Peridigm](https://github.com/peridigm/peridigm/blob/master/src/contact/Peridigm_ShortRangeForceContactModel.cpp)

"""
function compute_contact_model(cg, p::PenaltyContactParams, contact,
                               compute_master_force_density::Function,
                               compute_slave_force_density::Function)
    contact_dict = Data_Manager.get_contact_dict(cg)
    # function barrier on the dof: the forces below are static vectors
    with_dof(D -> penalty_forces!(contact_dict, p, contact, compute_master_force_density,
                                  compute_slave_force_density, D))
end

function penalty_forces!(contact_dict, p::PenaltyContactParams, contact,
                         compute_master_force_density::Function,
                         compute_slave_force_density::Function, ::Val{D}) where {D}
    contact_stiffness::Float64 = p.contact_stiffness
    contact_radius::Float64 = contact.base.contact_radius
    symmetry::String = contact.base.symmetry
    #Data_Manager.get_symmetry() # TODO store in materials the information
    for (master_id, contact) in pairs(contact_dict)
        distances = contact["Distances"]::Vector{Float64}
        normals = contact["Normals"]::Matrix{Float64}
        slaves = contact["Slaves"]::Vector{Int64}
        for id in 1:(contact["nSlaves"]::Int64)
            slave_id = slaves[id]
            horizon::Float64 = get_shared_horizon(slave_id) # needed to get the correct contact horizon
            # TODO symmetry needed
            if symmetry == "plane stress"
                stiffness = 9 / (pi * horizon^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
            elseif symmetry == "plane strain"
                stiffness = 48 / (5 * pi * horizon^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
            else
                stiffness = 9 / (pi * horizon^5)  # -> from Peridigm
                #stiffness = 12 / (pi * horizon^4)  # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
            end

            distance::Float64 = distances[id]
            normal = SVector{D,Float64}(ntuple(d -> normals[id, d], Val(D)))

            temp = contact_stiffness * stiffness *
                   (contact_radius - distance)
            normal_force = temp .* normal
            friction_id,
            friction_slave_id = compute_friction(id, slave_id,
                                                 p.friction_coefficient,
                                                 normal_force, normal)

            compute_master_force_density(master_id,
                                         slave_id,
                                         normal_force .+
                                         friction_id)
            compute_slave_force_density(slave_id,
                                        master_id,
                                        normal_force .+
                                        friction_slave_id)
        end
    end
end
"""

code was taken from Peridigm

"""
function compute_friction(id, slave_id, friction_coefficient,
                          normal_force::SVector{D,Float64},
                          normal::SVector{D,Float64}) where {D}
    no_friction = zero(SVector{D,Float64})
    if friction_coefficient == 0
        return no_friction, no_friction
    end

    mapping = Data_Manager.get_exchange_id_to_local_id()

    if isnothing(get(mapping, id, nothing)) || isnothing(get(mapping, slave_id, nothing))
        return no_friction, no_friction
    end

    velocity::NodeVectorField{Float64} = Data_Manager.get_field("Velocity", "NP1")
    velo_master = SVector{D,Float64}(ntuple(d -> velocity[mapping[id], d], Val(D)))
    velo_slave = SVector{D,Float64}(ntuple(d -> velocity[mapping[slave_id], d], Val(D)))

    current_dot_normal = dot(velo_master, normal)
    current_dot_neighbor = dot(velo_slave, normal)
    velo_id = velo_master .- current_dot_normal .* normal
    velo_slave_id = velo_slave .-
                    current_dot_neighbor .* normal
    vel_cm = 0.5 .* (velo_slave_id + velo_id)
    velo_id = velo_id .- vel_cm
    velo_slave_id = velo_slave_id .- vel_cm

    norm_vel_id = norm(velo_id)
    norm_velo_slave_id = norm(velo_slave_id)

    friction_id = norm_vel_id != 0 ?
                  -friction_coefficient * norm(normal_force) .* velo_id ./ norm_vel_id :
                  no_friction
    friction_slave_id = norm_velo_slave_id != 0 ?
                        -friction_coefficient * norm(normal_force) .* velo_slave_id ./
                        norm_velo_slave_id : no_friction
    return friction_id, friction_slave_id
end
end
