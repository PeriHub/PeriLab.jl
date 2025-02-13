# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using CSV
using DataFrames

function get_header(filename::Union{String,AbstractString})
    file = open(filename, "r")
    header_line = 0
    for line in eachline(file)#
        header_line += 1
        if contains(line, "header:")
            close(file)
            return header_line, convert(Vector{String}, split(line)[2:end])
        end
    end
    @error "No header exists in $filename. Please insert 'header: global_id' above the first node"
    return
end

function main(file, nodesets)
    header_line, header = get_header(file)
    global mesh = CSV.read(
        file,
        DataFrame;
        delim = " ",
        ignorerepeated = true,
        header = header,
        skipto = header_line + 1,
        comment = "#",
    )
    for nodeset in nodesets
        node_set_file = open(nodeset["file"], "w")
        println(node_set_file, "header: global_id")

        for id = 1:size(mesh, 1)
            global i = id
            if eval(Meta.parse(nodeset["function"]))
                println(node_set_file, "$id")
            end
        end
        close(node_set_file)
    end


end

file = "Zugstab.txt"
nodesets = [
    Dict("file" => "nodeset_1.txt", "function" => "mesh[!, \"x\"][i] < 30"),
    Dict("file" => "nodeset_2.txt", "function" => "mesh[!, \"x\"][i] > 220"),
]

main(file, nodesets)
