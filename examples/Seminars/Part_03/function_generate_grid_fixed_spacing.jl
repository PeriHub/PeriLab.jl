# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function generate_grid_fixed_spacing(filename, x_from, x_to, y_from, y_to, dx; block_id = 1)
    num_points_x = Int(floor((x_to - x_from) / dx)) + 1
    num_points_y = Int(floor((y_to - y_from) / dx)) + 1
    volume = dx * dx

    open(filename * ".txt", "w") do io
        println(io, "header: x y block_id volume number")
        for j in 0:(num_points_y - 1)
            y = y_from + j * dx
            for i in 0:(num_points_x - 1)
                x = x_from + i * dx
                println(io, "$(x) $(y) $(block_id) $(volume) $i")
            end
        end
    end
    open(filename * "_ns1.txt", "w") do io
        println(io, "header: global_id")
        for j in 0:(num_points_y - 1)
            y = 1 + j * num_points_x
            println(io, "$y")
        end
    end
    open(filename * "_ns2.txt", "w") do io
        println(io, "header: global_id")
        for j in 1:(num_points_y)
            y = j * num_points_x
            println(io, "$y")
        end
    end
end

generate_grid_fixed_spacing("seminar_03", 0, 0.05, 0, 0.01, 0.0025)
