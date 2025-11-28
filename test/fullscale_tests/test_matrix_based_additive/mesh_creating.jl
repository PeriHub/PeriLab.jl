# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Generate points for progressive activation (printing process)

Points are activated row by row from y=0 upwards, simulating
a layer-by-layer printing process.
"""

using Printf

# Parameters
x_min = 0.0
x_max = 0.01
y_min = 0.0
y_max = 0.005

dx = 0.0005  # Spacing in x-direction
dy = dx  # Spacing in y-direction (layer height)

volume = 2.5e-07  # Volume per point
block_id = 1  # Block ID for all points

# Activation parameters
activation_speed = 50.0  # Speed in x-direction (units/time)
layer_time = 0.1        # Time to switch between layers

# Generate grid points
x_coords = x_min:dx:x_max
y_coords = y_min:dy:y_max

n_x = length(x_coords)
n_y = length(y_coords)

println("Generating $(n_x * n_y) points...")
println("Grid: $n_x × $n_y")
println()

# Open output file
output_file = "mesh.txt"
open(output_file, "w") do f
    # Write header
    println(f, "header: x y block_id volume Activation_Time")
    activation_time=0
    # Generate points layer by layer (y-direction)
    for (j_idx, y) in enumerate(y_coords)
        # Time to reach this layer

        # Points in this layer (x-direction)
        for (i_idx, x) in enumerate(x_coords)
            # Activation time: layer start + time to reach x position

            activation_time += 1 / activation_speed

            # Write point data
            @printf(f, "%.15f %.15f %d %.15f %.15f\n",
                    x, y, block_id, volume, activation_time)
        end
    end
end

println("Points written to: $output_file")
println("Total points: $(n_x * n_y)")
println("Total activation time: $(@sprintf("%.3f", (n_y-1)*layer_time + (x_max-x_min)/activation_speed))")
