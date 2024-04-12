# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus
using Plots, Polynomials
using CSV, DataFrames

function main()
    if length(ARGS) < 1
        println("Usage: julia script.jl <filename>")
        return
    end

    # Extract filename from command-line arguments
    file_name = ARGS[1]

    exo = ExodusDatabase(file_name, "r")

    coords = read_coordinates(exo)
    damage = read_values(exo, NodalVariable, 12, "Damage")

    # Find indices of damaged coordinates
    damaged_indices = findall(x -> x > 0.3, damage)

    # Extract damaged coordinates
    damaged_coords = coords[:, damaged_indices]

    # Extract x and y values
    x_values = damaged_coords[1, :]
    y_values = damaged_coords[2, :]

    origin = (25, 0.0)
    min_ind = argmin(x_values)
    offset = (x_values[min_ind] - origin[1], y_values[min_ind] - origin[2])
    @info offset

    x_values .-= offset[1]
    # y_values .-= offset[2]

    times = read_times(exo)

    # f1 = fit(x_values, y_values, 2; weights=damage[damaged_indices])
    # f2 = fit(x_values, y_values, 4; weights=damage[damaged_indices])
    f3 = fit(x_values, y_values, 8; weights=damage[damaged_indices])
    # f4 = fit(x_values, y_values, 12; weights=damage[damaged_indices])
    # f5 = fit(x_values, y_values, 24; weights=damage[damaged_indices])

    scatter(x_values, y_values, zcolor=damage[damaged_indices], label="Damage")

    # plot!(f1, extrema(x_values)..., label="2")
    # plot!(f2, extrema(x_values)..., label="4")
    plot!(f3, extrema(x_values)..., label="8")
    # plot!(f4, extrema(x_values)..., label="12")
    # plot!(f5, extrema(x_values)..., label="24")
    savefig("plot.png")

    @info f3

    # Create a DataFrame to store the information
    range = vcat(minimum(x_values):0.1:maximum(x_values))
    df = DataFrame(
        x=range,
        f3=f3.(range),
    )
    CSV.write("plot.csv", df)
end

main()