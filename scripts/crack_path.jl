# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus
using Plots, Polynomials
using CSV, DataFrame

function main()
    if length(ARGS) < 1
        println("Usage: julia script.jl <filename>")
        return
    end

    # Extract filename from command-line arguments
    file_name = ARGS[1]

    exo = ExodusDatabase(file_name, "r")

    coords = read_coordinates(exo)
    times = read_times(exo)
    damage = read_values(exo, NodalVariable, length(times), "Damage")

    # Find indices of damaged coordinates
    damaged_indices = findall(x -> x > 0.2, damage)

    # Extract damaged coordinates
    damaged_coords = coords[:, damaged_indices]

    # Extract x and y values
    x_values = damaged_coords[1, :]
    y_values = damaged_coords[2, :]

    f1 = fit(x_values, y_values, 2; weights=damage[damaged_indices])
    f2 = fit(x_values, y_values, 4; weights=damage[damaged_indices])
    f3 = fit(x_values, y_values, 8; weights=damage[damaged_indices])
    f4 = fit(x_values, y_values, 12; weights=damage[damaged_indices])

    scatter(x_values, y_values, zcolor=damage[damaged_indices], label="Damage")

    plot!(f1, extrema(x_values)..., label="2")
    plot!(f2, extrema(x_values)..., label="4")
    plot!(f3, extrema(x_values)..., label="8")
    plot!(f4, extrema(x_values)..., label="12")
    savefig("plot.png")

    # Create a DataFrame to store the information
    range = vcat(minimum(x_values):1:maximum(x_values))
    df = DataFrame(
        x=range,
        f1=f1.(range),
    )
    CSV.write("plot.csv", df)
end

main()