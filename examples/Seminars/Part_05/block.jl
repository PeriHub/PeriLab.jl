# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function generate_grid_fixed_spacing(filename, x_from, x_to, y_from, y_to, dx; block_id = 1)
    num_points_x = Int(floor((x_to - x_from) / dx)) + 1
    num_points_y = Int(floor((y_to - y_from) / dx)) + 1
    volume = dx * dx
    top_list = []
    count = 0
    open(filename * ".txt", "w") do io
        println(io, "header: x y block_id volume number")
        for j in 0:(num_points_y - 1)
            y = y_from + j * dx
            if iseven(j)
                for i in 0:(num_points_x - 2)
                    x = x_from + i * dx + dx / 2
                    count += 1
                    println(io, "$(x) $(y) $(block_id) $(volume) $i")
                    if j == (num_points_y - 1)
                        push!(top_list, count)
                    end
                end
            else
                for i in 0:(num_points_x - 1)
                    x = x_from + i * dx
                    count += 1
                    println(io, "$(x) $(y) $(block_id) $(volume) $i")
                    if j == (num_points_y - 1)
                        push!(top_list, count)
                    end
                end
            end
        end
    end
    open(filename * "_ns1.txt", "w") do io
        println(io, "header: global_id")
        for j in 1:num_points_x
            println(io, "$j")
        end
    end
    open(filename * "_ns2.txt", "w") do io
        println(io, "header: global_id")
        for j in top_list
            println(io, "$j")
        end
    end
end

generate_grid_fixed_spacing("seminar_05", 0, 0.01, 0, 0.01, 0.0005)
