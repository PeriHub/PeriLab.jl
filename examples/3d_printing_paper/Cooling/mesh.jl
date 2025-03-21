# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Function to generate a simple mesh file with node sets
function write_mesh(x_end, y_end, dx, dy)

    # Create or open the file for writing
    filename = "mesh_000025.txt"
    mesh_file = open(filename, "w")

    # Calculate volume
    volume = dx * dy

    # Open node set files
    # node_set_file1 = open("BCall.txt", "w")
    # node_set_file2 = open("ns_simpleMesh_2.txt", "w")

    # Write header to the files
    println(mesh_file, "header: x y block_id volume")
    # println(node_set_file1, "header: global_id")
    # println(node_set_file2, "header: global_id")

    id = 1
    node_set = []

    # Loop through x and y coordinates to create the mesh
    for x in 0:dx:x_end
        for y in 0:dy:y_end
            block_id = 1
            # println(node_set_file1, "$id")

            # Determine block_id based on x coordinate
            # if x <= 2
            #     block_id = 2
            # elseif x >= 8
            #     block_id = 3
            #     println(node_set_file2, "$id")
            # end

            # Write mesh information to the file
            println(mesh_file, "$x $y $block_id $volume")
            id += 1
        end
    end

    # Close the files
    close(mesh_file)
    # close(node_set_file1)
    # close(node_set_file2)

    # Print success message
    println("File '$filename' created successfully.")
end

# Call the function with specified parameters
write_mesh(0.01, 0.01, 0.000125, 0.000125)
