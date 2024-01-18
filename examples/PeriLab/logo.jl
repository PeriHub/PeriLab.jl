# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Images

function write_mesh()
    # Create or open the file for writing
    filename = "PeriLab.txt"
    file = open(filename, "w")
    node_set_file1 = open("ns_PeriLab_1.txt", "w")
    node_set_file2 = open("ns_PeriLab_2.txt", "w")

    # Write the header to the file
    println(file, "header: x y block_id volume")
    println(node_set_file1, "header: global_id")
    println(node_set_file2, "header: global_id")

    # Text to be converted to pixels
    # text = "PeriLab"

    # Convert text to image
    img = load("logo.png")
    @info size(img)
    id = 1
    d_x = 10
    d_y = 10
    scale = 0.01
    node_set = []
    # Write the coordinates and parameters to the file
    for x in 1:d_x:size(img, 2)
        for y in 1:d_y:size(img, 1)
            # @info img[y, x]
            volume = d_x * scale * d_y * scale
            y_inv = (size(img, 1) - y) * scale
            x_new = x * scale
            if x_new < 0.5
                println(node_set_file1, "$id")
            elseif x_new > 17.34
                println(node_set_file2, "$id")
            end
            if img[y, x] != RGBA{N0f8}(1.0, 1.0, 1.0, 1.0)  # Check if pixel is part of the text
                block_id = 1
                println(file, "$x_new $y_inv $block_id $volume")
            elseif x_new <= 0.7
                block_id = 3
                println(file, "$x_new $y_inv $block_id $volume")
            elseif x_new >= 17
                block_id = 4
                println(file, "$x_new $y_inv $block_id $volume")
            else
                block_id = 2
                println(file, "$x_new $y_inv $block_id $volume")
            end
            id += 1
        end
    end

    # Close the file
    close(file)
    close(node_set_file1)
    close(node_set_file2)

    println("File '$filename' created successfully.")
end

write_mesh()