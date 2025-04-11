# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Penalty_model

function init_contact_model(datamanager, params)
end

function compute_contact_model(datamanager, params, distances, projection, contact, forces,
                               main_id, surface_subset)

    # distances are given by the global ids -> projection must be made outside
    # distances, projection,contact, forces provided in the correct style

    # distance, projected distance, norm distance, main_id (master or slave), surface_subset,

    # distribution between master and slave is made bei params -> must be made, because they can be on different cores
    # neighbor volume might not exist here (slave  or master)
    # ids are numbered lists

    # for (id, point_id) in enumerate(cont)
    # if master
    # force[point_id] += dist[i]*stiffness
    # if slave
    # force[point_id] += dist[i]*stiffness
    #force[subset_id] += dist[i]*stiffness
    # end

end

end
