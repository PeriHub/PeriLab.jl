# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Images

function write_mesh(filename, dx, dy, scale)

    img = load(filename)
    @info "Image size: " size(img)
    #get unique colors
    colors = unique(img)
    #remove white
    colors = colors[colors.!=RGBA(1.0, 1.0, 1.0, 1.0)]
    num_blocks = length(colors)
    @info "Found $(num_blocks) block(s)"

    # Create or open the file for writing
    name = split(filename, ".")[1]
    file = open(name * ".txt", "w")

    node_set_files = []
    for id in 1:num_blocks
        push!(node_set_files, open("ns_" * name * "_" * string(id) * ".txt", "w"))
        println(node_set_files[id], "header: global_id")
        id += 1
    end

    # Write the header to the file
    println(file, "header: x y block_id volume")
    id = 1
    # Write the coordinates and parameters to the file
    for x in 1:dx:size(img, 2)
        for y in size(img, 1):-dy:1
            # @info img[y, x]
            volume = dx * scale * dy * scale
            y_inv = (size(img, 1) - y) * scale
            x_new = x * scale
            # get index of color in colors array
            block_id = findfirst(color -> color == img[y, x], colors)
            if !isnothing(block_id)
                println(node_set_files[block_id], "$id")
                println(file, "$x_new $y_inv $block_id $volume")
                id += 1
            end
        end
    end
    @info "Found $(id - 1) node(s)"
    # Close the file
    close(file)
    for id in 1:num_blocks
        close(node_set_files[id])
    end

    @info "Model '$name' created successfully."
end

filename = "Image.png"
dx = 10
dy = 10
scale = 0.01
write_mesh(filename, dx, dy, scale)