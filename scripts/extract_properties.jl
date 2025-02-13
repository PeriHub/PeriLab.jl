# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function main()
    file = open("input.inp", "r")
    lines = readlines(file)
    close(file)
    properties = []
    read = false
    for line in lines
        if occursin("*User Material", line)
            read = true
            continue
        end
        if read
            if occursin("*", line)
                break
            end
            id = 1
            values = split(line, ",")
            for value in values
                push!(properties, parse(Float64, strip(value)))
                id += 1
            end
            for i = 1:8-id
                push!(properties, nothing)
            end
        end
    end
    @info "Found $(length(properties)) properties"
    for (id, prop) in enumerate(properties)
        if isnothing(prop)
            continue
        end
        @info "Property_$(id): $(prop)"
    end
end

main()
