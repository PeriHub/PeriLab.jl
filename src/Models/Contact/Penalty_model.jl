# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Penalty_model

include("../../Support/Helpers.jl")
import .Helpers: get_shared_horizon

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
    return datamanager
end

function compute_contact_model(datamanager, cm, params, compute_master_force_density,
                               compute_slave_force_density)
    contact_dict = datamanager.get_contact_dict(cm)
    contact_stiffness = params["Contact Stiffness"]
    contact_radius = params["Contact Radius"]
    #datamanager.get_symmetry() # TODO store in materials the information
    for (master_id, contact) in pairs(contact_dict)
        for (id, slave_id) in enumerate(contact["Slaves"])
            horizon = get_shared_horizon(datamanager, slave_id) # needed to get the correct contact horizon
            # TODO symmetry needed
            stiffness = 48.0 / 5.0 * contact_stiffness / (pi * horizon^3)

            @views distance = contact["Distances"][id]
            temp = (contact_radius - distance) / horizon

            @views normal = contact["Normals"][id]
            compute_master_force_density(datamanager, master_id, slave_id,
                                         stiffness * temp .*
                                         normal)
            compute_slave_force_density(datamanager, slave_id, master_id,
                                        stiffness * temp .*
                                        normal)
        end
    end
    return datamanager
end

end
