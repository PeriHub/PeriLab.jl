# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function generate_grid_fixed_spacing(filename, x_from, x_to, y_from, y_to, dx; block_id = 1)
    num_points_x = Int(floor((x_to - x_from) / dx)) + 1
    num_points_y = Int(floor((y_to - y_from) / dx)) + 1
    volume = dx * dx

    open(filename, "w") do io
        println(io, "header: x y block_id volume number")
        for j in 0:(num_points_y - 1)
            y = y_from + j * dx
            for i in 0:(num_points_x - 1)
                x = x_from + i * dx
                println(io, "$(x) $(y) $(block_id) $(volume) $i")
            end
        end
    end
end

generate_grid_fixed_spacing("seminar_02.txt", 0, 1, -1, 3, 0.2)
