# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using CSV, DataFrames
using LazyGrids
using ProgressBars
using LinearAlgebra
using Plots

function closest_point_to_vector(start_point::Vector{Float64}, end_point::Vector{Float64}, point::Vector{Float64})
    # Calculate the direction vector of the line segment
    dir = normalize(end_point - start_point)

    # Calculate the distance from the point to the line segment
    distance_along_line = dot(point - start_point, dir) / dot(dir, dir)
    # distance_along_line = clamp(distance_along_line, 0.0, 1.0)

    # Calculate the closest point on the line segment
    closest_point = start_point + (dir * distance_along_line)

    # Calculate the distance from the closest point to the point
    distance_to_closest_point = sqrt((point[1] - closest_point[1])^2 + (point[2] - closest_point[2])^2)

    return distance_along_line, distance_to_closest_point
end

function write_mesh(event_file, dx, dt, scale, width)

    events = CSV.read(event_file, DataFrame; delim=",", comment="#")

    time_event = events[:, 1]
    x_event = events[:, 2]
    y_event = events[:, 3]
    z_event = events[:, 4]
    layers = unique(z_event)
    extrude_event = events[:, 5]

    height = minimum(layers)

    dx_value = [dx * scale, dx * scale, height / 2]

    volume = dx_value[1] * dx_value[2] * dx_value[3]

    x_min = minimum(x_event) - width / 2
    x_max = maximum(x_event) + width / 2

    y_min = minimum(y_event) - width / 2
    y_max = maximum(y_event) + width / 2

    (xg, yg) = ndgrid(x_min:dx_value[1]:x_max, y_min:dx_value[2]:y_max)

    grid_x = reduce(vcat, xg)
    grid_y = reduce(vcat, yg)

    previous_time = 0
    previous_x = 0
    previous_y = 0
    previous_z = 0
    previous_extruding = 0
    x_peri = []
    y_peri = []
    used_ids = []
    p = Plots.plot()

    mesh_df = DataFrame(x=Float64[], y=Float64[], z=Float64[], k=Int64[], volume=Float64[], time=Float64[])

    iter = ProgressBar(eachindex(time_event))
    for i in iter

        time = time_event[i]
        x = x_event[i]
        y = y_event[i]
        z = z_event[i]
        extruding = extrude_event[i]

        if previous_extruding == 1
            distance = norm([x, y] - [previous_x, previous_y])
            v = distance / (time - previous_time)

            p = plot!(p, [previous_x, x], [previous_y, y], label="Plot", lc=:black, lw=1)

            for j in eachindex(grid_x)
                px = grid_x[j]
                py = grid_y[j]
                if !(j in used_ids)
                    distance_along_line, distance_to_closest_point = closest_point_to_vector([previous_x, previous_y], [x, y], [px, py])
                    if distance_to_closest_point <= width / 2 && distance_along_line <= distance && distance_along_line >= 0
                        time_to_activation = distance_along_line / v
                        push!(used_ids, j)
                        push!(mesh_df, [px, py, z, 1, volume, time_to_activation + previous_time])
                        push!(x_peri, px)
                        push!(y_peri, py)
                    end
                end
            end
        end

        if i == length(time_event) || z != z_event[i+1]
            p = scatter!(p, x_peri, y_peri, title="Layer" * string(z), xlabel="X", ylabel="Y", ma=0.5, ms=1)
            savefig(p, "Output/layer" * string(z) * ".svg")
            x_peri = []
            y_peri = []
            used_ids = []
            p = Plots.plot()
            extruding = 0
        end

        previous_x = x
        previous_y = y
        previous_z = z
        previous_time = time
        previous_extruding = extruding

    end

    txt_file = joinpath("Output", replace(event_file, ".inp" => ".txt"))
    write(txt_file, "header: x y z block_id volume Activation_Time\n")
    CSV.write(txt_file, mesh_df; delim=' ', append=true)
end

event_file = "L-Angle_v2_EventSeries.inp"
dx = 0.2
dt = 0.002
scale = 1
width = 0.5

write_mesh(event_file, dx, dt, scale, width)