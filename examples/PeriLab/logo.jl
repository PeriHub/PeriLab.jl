# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Images

function write_mesh()
    # Create or open the file for writing
    filename = "PeriLab.txt"
    file = open(filename, "w")
    node_set_file = open("ns_PeriLab.txt", "w")

    # Write the header to the file
    println(file, "header: x y block_id volume")
    println(node_set_file, "header: global_id")

    # Text to be converted to pixels
    # text = "PeriLab"

    # Convert text to image
    img = load("logo.png")
    @info size(img)
    id = 1
    d_x = 10
    d_y = 10
    node_set = []
    # Write the coordinates and parameters to the file
    for x in 1:d_x:size(img, 2)
        for y in 1:d_y:size(img, 1)
            # @info img[y, x]
            if img[y, x] != RGBA{N0f8}(1.0, 1.0, 1.0, 1.0)  # Check if pixel is part of the text
                block_id = 1
                volume = d_x * d_y
                y_inv = size(img, 1) - y
                if x < 1150 && x > 1000 && y_inv > 400
                    block_id = 2
                    println(node_set_file, "$id")
                end
                println(file, "$x $y_inv $block_id $volume")
                id += 1
            end
        end
    end

    # Close the file
    close(file)
    close(node_set_file)

    println("File '$filename' created successfully.")
end

write_mesh()