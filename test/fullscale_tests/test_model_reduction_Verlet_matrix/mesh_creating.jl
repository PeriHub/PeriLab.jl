# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Printf
"""
Generate points for progressive activation (printing process)
Points are activated row by row from y=0 upwards, simulating
a layer-by-layer printing process with multiple blocks in x-direction.
"""
# Parameters
x_min = 0.0
x_max = 0.01
y_min = 0.0
y_max = 0.002
z_min = 0.0
z_max = 0.002
dx = 0.0004  # Spacing in x-direction
dy = dx      # Spacing in y-direction (layer height)
dz = dx      # Spacing in z-direction
volume = 2.5e-07  # Volume per point

# Activation parameters
activation_speed = 50.0  # Speed in x-direction (units/time)
layer_time = 0.1         # Time to switch between layers

# Block definition in x-direction
# Define block boundaries and their IDs
block_definitions = [
    (x_start = x_min, x_end = 0.002, block_id = 1),
    (x_start = 0.002, x_end = 0.004, block_id = 2),
    (x_start = 0.004, x_end = 0.006, block_id = 1),
    (x_start = 0.006, x_end = 0.007, block_id = 2),
    (x_start = 0.007, x_end = 0.009, block_id = 3),
    (x_start = 0.009, x_end = x_max, block_id = 1)
]

"""
	get_block_id(x, block_defs)

Determine block ID based on x-coordinate.
"""
function get_block_id(x, block_defs)
    for block in block_defs
        if block.x_start <= x < block.x_end ||
           (abs(x - block.x_end) < 1e-10 && block == last(block_defs))
            return block.block_id
        end
    end
    return 1  # Default block
end

"""
	is_in_first_block(x, block_defs)

Check if point is in first block.
"""
function is_in_first_block(x, block_defs)
    first_block = first(block_defs)
    return first_block.x_start <= x <= first_block.x_end
end

"""
	is_in_last_block(x, block_defs)

Check if point is in last block.
"""
function is_in_last_block(x, block_defs)
    last_block = last(block_defs)
    return last_block.x_start <= x <= last_block.x_end
end

# Generate grid points
x_coords = x_min:dx:x_max
y_coords = y_min:dy:y_max
z_coords = z_min:dz:z_max

# ============================================
# 2D MESH GENERATION
# ============================================
n_x = length(x_coords)
n_y = length(y_coords)
println("=" ^ 50)
println("2D MESH GENERATION")
println("=" ^ 50)
println("Generating $(n_x * n_y) points...")
println("Grid: $n_x × $n_y")
println("Blocks in x-direction: $(length(block_definitions))")
println()

output_file = "mesh_2d.txt"
nset_1_file = "nset_1_2d.txt"
nset_2_file = "nset_2_2d.txt"

# Open all files
open(output_file, "w") do f
    open(nset_1_file, "w") do f1
        open(nset_2_file, "w") do f2
            # Write headers
            println(f, "header: x y block_id volume")
            println(f1, "header: global_id")
            println(f2, "header: global_id")

            node_id = 1

            # Generate points layer by layer (y-direction)
            for (j_idx, y) in enumerate(y_coords)
                # Points in this layer (x-direction)
                for (i_idx, x) in enumerate(x_coords)
                    # Determine block ID based on x position
                    block_id = get_block_id(x, block_definitions)

                    # Write point data to main file
                    @printf(f, "%.15f %.15f %d %.15e\n",
                            x, y, block_id, volume)

                    # Write to node set files if in first or last block
                    if is_in_first_block(x, block_definitions)
                        @printf(f1, "%d\n", node_id)
                    end

                    if is_in_last_block(x, block_definitions)
                        @printf(f2, "%d\n", node_id)
                    end

                    node_id += 1
                end
            end
        end
    end
end

println("2D Points written to: $output_file")
println("Node set 1 (first block) written to: $nset_1_file")
println("Node set 2 (last block) written to: $nset_2_file")
println("Total points: $(n_x * n_y)")
println()

# ============================================
# 3D MESH GENERATION
# ============================================
n_z = length(z_coords)
n_points_3d = n_x * n_y * n_z
println("=" ^ 50)
println("3D MESH GENERATION")
println("=" ^ 50)
println("Generating $n_points_3d points...")
println("Grid: $n_x × $n_y × $n_z")
println("Blocks in x-direction: $(length(block_definitions))")
println()

output_file = "mesh_3d.txt"
nset_1_file = "nset_1_3d.txt"
nset_2_file = "nset_2_3d.txt"

open(output_file, "w") do f
    open(nset_1_file, "w") do f1
        open(nset_2_file, "w") do f2
            # Write headers
            println(f, "header: x y z block_id volume")
            println(f1, "header: global_id")
            println(f2, "header: global_id")

            node_id = 1

            # Generate points layer by layer (z-direction, then y-direction)
            for (k_idx, z) in enumerate(z_coords)
                for (j_idx, y) in enumerate(y_coords)
                    # Points in this layer (x-direction)
                    for (i_idx, x) in enumerate(x_coords)
                        # Determine block ID based on x position
                        block_id = get_block_id(x, block_definitions)

                        # Write point data to main file
                        @printf(f, "%.15f %.15f %.15f %d %.15e\n",
                                x, y, z, block_id, volume)

                        # Write to node set files if in first or last block
                        if is_in_first_block(x, block_definitions)
                            @printf(f1, "%d\n", node_id)
                        end

                        if is_in_last_block(x, block_definitions)
                            @printf(f2, "%d\n", node_id)
                        end

                        node_id += 1
                    end
                end
            end
        end
    end
end

println("3D Points written to: $output_file")
println("Node set 1 (first block) written to: $nset_1_file")
println("Node set 2 (last block) written to: $nset_2_file")
println("Total points: $n_points_3d")
println()

# ============================================
# SUMMARY
# ============================================
println("=" ^ 50)
println("SUMMARY")
println("=" ^ 50)
println("Block definitions:")
for (i, block) in enumerate(block_definitions)
    println("  Block $i: x ∈ [$(block.x_start), $(block.x_end)] → ID=$(block.block_id)")
end
println()
println("First block: x ∈ [$(first(block_definitions).x_start), $(first(block_definitions).x_end)]")
println("Last block:  x ∈ [$(last(block_definitions).x_start), $(last(block_definitions).x_end)]")
