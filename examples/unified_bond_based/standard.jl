# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_node_header()
    return "header: global_id\n"
end
function get_mesh_header(dof = 2)
    if dof == 2
        return "header: x y block_id volume\n"
    else
        return "header: x y z block_id volume\n"
    end
end
