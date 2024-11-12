# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using CSV, DataFrames
using LazyGrids
using ProgressBars
using LinearAlgebra
using Plots

function fastdot(a, b)
    c = 0.0
    @inbounds @simd for i ∈ eachindex(a, b)
        c += a[i] * b[i]
    end
    c
end

function sub_in_place!(C::Vector{T}, A::Vector{T}, B::Vector{T}) where {T<:Number}
    @assert length(C) == length(A) == length(B)

    @inbounds for i ∈ eachindex(A)
        C[i] = A[i] - B[i]
    end
end

function normalize_in_place!(B::Vector{T}, A::Vector{T}) where {T<:Number}
    nrm = norm(A)
    @inbounds for i ∈ eachindex(A)
        B[i] = A[i] / nrm
    end
end

function closest_point_to_vector(
    start_point::Vector{Float64},
    dir::Vector{Float64},
    point::Vector{Float64},
    closest_point::Vector{Float64},
    point_diff::Vector{Float64},
)

    # Calculate the distance from the point to the line segment
    distance_along_line = dot(point_diff, dir) / dot(dir, dir)

    # Calculate the closest point on the line segment
    closest_point .= start_point .+ (dir .* distance_along_line)

    # Calculate the distance from the closest point to the point
    distance_to_closest_point =
        sqrt((point[1] - closest_point[1])^2 + (point[2] - closest_point[2])^2)

    return distance_along_line, distance_to_closest_point
end

function write_mesh(event_file, dx, dt, scale, width, plot, plot_all = false)
    if !isdir("Output")
        mkdir("Output")
    end

    events = CSV.read(event_file, DataFrame; delim = ",", comment = "#")

    time_event = events[:, 1]
    x_event = events[:, 2]
    y_event = events[:, 3]
    z_event = events[:, 4]
    layers = unique(z_event)
    extrude_event = events[:, 5]

    height = minimum(layers)

    dx_value = [dx * scale, dx * scale, (height / 2) * scale]

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
    used_ids = zeros(length(grid_x))
    p = Plots.plot()

    mesh_df = DataFrame(
        x = Float64[],
        y = Float64[],
        z = Float64[],
        k = Int64[],
        volume = Float64[],
        time = Float64[],
    )
    closest_point = zeros(2)
    dir = zeros(2)
    start_point = zeros(2)
    point = zeros(2)
    point_diff = zeros(2)

    iter = ProgressBar(eachindex(time_event))
    for i in iter

        time = time_event[i]
        x = x_event[i]
        y = y_event[i]
        z = z_event[i]
        extruding = extrude_event[i]

        if previous_extruding == 1
            start_point[1] = previous_x
            start_point[2] = previous_y
            point[1] = x
            point[2] = y
            sub_in_place!(point_diff, point, start_point)
            distance = norm(point_diff)
            v = distance / (time - previous_time)

            if plot
                p = plot!(
                    p,
                    [previous_x, x],
                    [previous_y, y],
                    # label = "Plot",
                    legend = false,
                    lc = :black,
                    lw = 1,
                )
            end

            normalize_in_place!(dir, point_diff)
            for j in eachindex(grid_x)
                point[1] = grid_x[j]
                point[2] = grid_y[j]
                sub_in_place!(point_diff, point, start_point)
                if used_ids[j] == 0
                    distance_along_line, distance_to_closest_point =
                        closest_point_to_vector(
                            start_point,
                            dir,
                            point,
                            closest_point,
                            point_diff,
                        )
                    if distance_to_closest_point <= width / 2 &&
                       distance_along_line <= distance &&
                       distance_along_line >= 0
                        time_to_activation = distance_along_line / v
                        used_ids[j] = 1
                        push!(
                            mesh_df,
                            [
                                point[1],
                                point[2],
                                z,
                                1,
                                volume,
                                time_to_activation + previous_time,
                            ],
                        )
                        if plot
                            push!(x_peri, point[1])
                            push!(y_peri, point[2])
                        end
                    end
                end
            end
        end

        if i == length(time_event) || z != z_event[i+1]
            if plot
                p = scatter!(
                    p,
                    x_peri,
                    y_peri,
                    title = "Layer" * string(z),
                    xlabel = "X",
                    ylabel = "Y",
                    ma = 0.5,
                    ms = 1,
                )
                savefig(p, "Output/layer" * string(z) * ".svg")
                x_peri = []
                y_peri = []
                p = Plots.plot()
                if !plot_all
                    plot = false
                end
            end
            fill!(used_ids, 0)
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
    CSV.write(txt_file, mesh_df; delim = ' ', append = true)
end

event_file = "L-Angle_v2_EventSeries.inp"
dx = 0.2
dt = 0.002
scale = 1
width = 0.5
plot = true

write_mesh(event_file, dx, dt, scale, width, plot)
