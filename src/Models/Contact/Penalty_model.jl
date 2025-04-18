# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Penalty_model

function init_contact_model(datamanager, params)
end

function compute_contact_model(datamanager, cm, params, compute_master_force_density,
                               compute_slave_force_density)
    contact_dict = datamanager.get_contact_dict()
    stiffness = params["Contact stiffness"]
    nlist = datamanager.get_nlist()
    contact_pairs = contact_dict["Pairs: Master-Slave"]
    distances = contact_dict["Pairs: Master-Slave"]

    for (contact_id, pair) in enumerate(contact_dict["Pairs: Master-Slave"])
        master_id = pair[1]
        distances = contact_dict["Distances"][contact_id]
        normals = contact_dict["Normals"][contact_id]

        for (pair_id, slave_id) in enumerate(pair[2])
            compute_master_force_density(datamanager, master_id,
                                         stiffness * distances * normals[pair_id, :])
            compute_slave_force_density(datamanager, slave_id,
                                        stiffness * distances * normals[pair_id, :])
        end
    end
    append!(contact_dict["Pairs: Master-Slave"], [(master_node, id)])
end

end
