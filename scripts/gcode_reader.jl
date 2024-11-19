# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
# Gcode functions taken from MIT Project GcodeParser.jl https://github.com/janvorisek/GcodeParser.jl

using Plots
using LinearAlgebra
using LazyGrids
using CSV, DataFrames
using ProgressBars
using NearestNeighbors
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--dx", "-x"
        help = "dx"
        arg_type = Float64
        default = 0.2
        "--dy", "-y"
        help = "dy"
        arg_type = Float64
        default = 0.2
        "--width", "-w"
        help = "width"
        arg_type = Float64
        default = 0.4
        "--scale", "-s"
        help = "scale"
        arg_type = Float64
        default = 1.0
        "--plot_enabled", "-p"
        help = "plot_enabled"
        arg_type = Bool
        default = true
        "filename"
        help = "filename"
        required = true
    end

    return parse_args(s)
end

function sub_in_place!(C::Vector{T}, A::Vector{T}, B::Vector{T}) where {T<:Number}
    @assert length(C) == length(A) == length(B)

    @inbounds for i ∈ eachindex(A)
        C[i] = A[i] - B[i]
    end
end

function normalize_in_place!(B::Vector{T}, A::Vector{T}) where {T<:Number}
    nrm = norm(A)
    @inbounds for i ∈ eachindex(A)
        B[i] = A[i] / nrm
    end
end
function closest_point_to_vector(
    start_point::Vector{Float64},
    dir::Vector{Float64},
    point::Vector{Float64},
    closest_point::Vector{Float64},
    point_diff::Vector{Float64},
)

    # Calculate the distance from the point to the line segment
    distance_along_line = dot(point_diff, dir) / dot(dir, dir)

    # Calculate the closest point on the line segment
    closest_point .= start_point .+ (dir .* distance_along_line)

    # Calculate the distance from the closest point to the point
    distance_to_closest_point =
        sqrt((point[1] - closest_point[1])^2 + (point[2] - closest_point[2])^2)

    return distance_along_line, distance_to_closest_point
end
"""
    stripComments(line::String)::String

Return a copy of string `line` with stripped comments inside parentheses and all characters after a semicolon.

This function also removes whitespace as it it not needed for further parsing.

# Examples
```julia-repl
julia> stripComments("G92 (G10(aaa)))) ((comment)G) Z0.2 ; this is a comment")
"G92Z0.2"
```
"""
function stripComments(line::String)::String
    re1 = r"\(.*\)"    # Remove anything inside the outer parentheses
    re2 = r"[^:]\;.*"  # Remove anything after a semi-colon to the end of the line, including preceding spaces

    line = replace(line, re1 => s"")
    line = replace(line, re2 => s"")
    line = filter(x -> !isspace(x), line) # Remove whitespace

    return line
end

"""
    parseLine(line::String, returnPair::Bool = true)::Array{Union{String,Pair{String,String}},1}

Parse a single line of g-code and return an array of `Pair{String,String}` or an array of `String` containing the parsed commands.

The first command usually defines what to do (ie. `G01` - linear interpolation) and following commands are the arguments (ie. `X 14.312`);

# Examples
```julia-repl
julia> parseLine("G10 X5.Y3. E6.")
4-element Array{Union{Pair{String,String}, String},1}:
 "G" => "10"
 "X" => "5."
 "Y" => "3."
 "E" => "6."
```

Return array of strings
```julia-repl
julia> parseLine("G10 X5.Y3. E6.", false)
4-element Array{Union{Pair{String,String}, String},1}:
 "G10"
 "X5."
 "Y3."
 "E6."
```
"""
function parseLine(
    line::String,
    returnPair::Bool = true,
)::Array{Union{String,Pair{String,String}},1}
    line = stripComments(line)

    # Match commands
    gcode_regex =
        r"/(%.*)|({.*)|((?:\$\$)|(?:\$[a-zA-Z0-9#]*))|([a-zA-Z][0-9\+\-\.]+)|(\*[0-9]+)/igm"

    # array of matched strings
    matches = collect(String(m.match) for m in eachmatch(gcode_regex, line))

    if returnPair
        return collect(first(m, 1) => last(m, length(m) - 1) for m in matches)
    end

    return matches
end

function parseFile(path::String, callbacks::Dict{String,Function}, dataObject)
    lines = readlines(path)
    # open(path) do f
    #     line = 1
    #     while !eof(f)
    #         x = readline(f);
    iter = ProgressBar(eachindex(lines))
    for i in iter
        x = lines[i]
        if occursin(";", x)
            if occursin(";TYPE:", x)
                command = split(x, ":")[2]
            elseif occursin(";Z:", x)
                command = "new_layer"
                z = parse(Float64, split(x, ":")[2])
                callbacks[command](z, dataObject)
                continue
            else
                command = x[2:end]
            end
            if haskey(callbacks, command)
                if dataObject === nothing
                    callbacks[command]()
                else
                    callbacks[command](dataObject)
                end
            end
        else
            cmds = parseLine(x)
            if length(cmds) == 0
                continue
            end

            letter = cmds[1].first
            number = cmds[1].second

            command = "$letter$number"
            if haskey(callbacks, command)
                if dataObject === nothing
                    callbacks[command](cmds)
                else
                    callbacks[command](cmds, dataObject)
                end
            end
        end
    end
    #         line += 1;
    #     end
    # end
end

function write_mesh(gcode_file, find_min_max, discretization, pd_mesh = Dict())

    # create any data object
    # it will be passed as a second parameter to your callbacks
    # here simple dictionary is used to store information during the print
    myPrinter = Dict{String,Any}()
    myPrinter["plot"] = Plots.plot()
    myPrinter["positioning"] = "absolute"
    myPrinter["x"] = 0.0
    myPrinter["y"] = 0.0
    myPrinter["previous_x"] = 0.0
    myPrinter["previous_y"] = 0.0
    myPrinter["z"] = 0.0
    myPrinter["e"] = 0.0
    myPrinter["filamentUsage"] = 0.0 # store total filament usage (printed length of filament)
    myPrinter["distanceMoved"] = 0.0 # store total distance moved
    myPrinter["time"] = 0.0
    myPrinter["previous_time"] = 0.0
    myPrinter["previous_extruding"] = false
    myPrinter["relevant_component"] = false
    myPrinter["find_min_max"] = find_min_max
    myPrinter["x_min"] = 1.e100
    myPrinter["x_max"] = 0.0
    myPrinter["y_min"] = 1.e100
    myPrinter["y_max"] = 0.0
    myPrinter["pd_mesh"] = pd_mesh
    myPrinter["layers"] = []

    # Setup a dictionary of callbacks for specified commands
    callbacks = Dict{String,Function}()
    callbacks["G0"] = move # just move the printhead
    callbacks["G1"] = extrude  # move the printhead as well as extrude material
    callbacks["Perimeter"] = switch_on
    callbacks["Solid infill"] = switch_on
    # callbacks["WIPE_START"] = switch_off
    # callbacks["WIPE_END"] = switch_on
    callbacks["Custom"] = switch_off
    callbacks["Skirt/Brim"] = switch_off
    callbacks["new_layer"] = new_layer

    # watch out for relative and absolute positioning
    callbacks["G90"] = (cmds, dataobject) -> dataobject["positioning"] = "absolute"
    callbacks["G91"] = (cmds, dataobject) -> dataobject["positioning"] = "relative"

    # parse g-code file and simulate print using our own callbacks and data object
    parseFile(gcode_file, callbacks, myPrinter)

    if find_min_max
        @info "Minimum and maximum values:"
        @info "X min/max:", myPrinter["x_min"], myPrinter["x_max"]
        @info "Y min/max:", myPrinter["y_min"], myPrinter["y_max"]
        layer_heights = round.(diff(myPrinter["layers"]); digits = 10)
        if !all(layer_heights .== layer_heights[1])
            @warn "Layer heights are not equal"
            @warn layer_heights
        end
        @info "Layers: $(myPrinter["layers"])"
        return ndgrid(
            myPrinter["x_min"]:discretization[1]:myPrinter["x_max"],
            myPrinter["y_min"]:discretization[2]:myPrinter["y_max"],
        ),
        layer_heights[1]
    end
    # Show printer data after print with some interesting stats
    # @show myPrinter;

    return
end

function move(cmds, dataobject)
    movement(cmds, dataobject)
    dataobject["previous_extruding"] = false
end

function check_min_max(dataobject, str)
    if dataobject[str] > dataobject[str*"_max"]
        dataobject[str*"_max"] = dataobject[str]
    end
    if dataobject[str] < dataobject[str*"_min"]
        dataobject[str*"_min"] = dataobject[str]
    end
end

"""
    movement(cmds, dataobject)

Example movement callback for `G0` and `G1` which calculates the total distance moved in all axes.

It is calculated by watching the `X`, `Y` and `Z` axes movement.
"""
function movement(cmds, dataobject)
    dataobject["previous_x"] = dataobject["x"]
    dataobject["previous_y"] = dataobject["y"]

    dx = 0.0
    dy = 0.0
    dz = 0.0

    x = findfirst((x -> lowercase(x.first) == "x"), cmds)
    if x !== nothing
        val = parse(Float64, cmds[x].second)

        if dataobject["positioning"] === "absolute"
            dx = val - dataobject["x"]
            dataobject["x"] = val
        else
            dx = val
        end
    end

    y = findfirst((x -> lowercase(x.first) == "y"), cmds)
    if y !== nothing
        val = parse(Float64, cmds[y].second)

        if dataobject["positioning"] === "absolute"
            dy = val - dataobject["y"]
            dataobject["y"] = val
        else
            dy = val
        end
    end

    z = findfirst((x -> lowercase(x.first) == "z"), cmds)
    if z !== nothing
        val = parse(Float64, cmds[z].second)

        if dataobject["positioning"] === "absolute"
            dz = val - dataobject["z"]
            dataobject["z"] = val
        else
            dz = val
        end
    end

    f = findfirst((x -> lowercase(x.first) == "f"), cmds)
    if f !== nothing
        val = parse(Float64, cmds[f].second)

        if dataobject["positioning"] === "absolute"
            dataobject["f"] = val
        end
    end
    distance = sqrt(dx * dx + dy * dy + dz * dz)
    dataobject["distanceMoved"] += distance
    dataobject["previous_time"] = dataobject["time"]
    dataobject["time"] += distance / dataobject["f"] * 60
end

"""
    extrude(cmds, dataobject)

Example extrusion callback for `G1` which calculates total length of filament extruded.

The extruded filament length is obtained by watching the `E` axis movement in the g-code file.
"""
function extrude(cmds, dataobject)
    movement(cmds, dataobject)

    # calculate used filament length
    e = findfirst((x -> lowercase(x.first) == "e"), cmds)
    if e !== nothing
        # Current E axis value
        e = parse(Float64, cmds[e].second)

        # Printed length of a current move
        if dataobject["positioning"] === "absolute"
            de = e - dataobject["e"]
            dataobject["e"] = e
        else
            de = e
        end

        if e <= 0.0
            return
        end

        # Used filament
        dataobject["filamentUsage"] += de
        # println(dataobject["filamentUsage"]);5
        if dataobject["previous_extruding"] && dataobject["relevant_component"]
            check_min_max(dataobject, "x")
            check_min_max(dataobject, "y")
            if !dataobject["find_min_max"]
                if dataobject["pd_mesh"]["plot_enabled"]
                    dataobject["plot"] = plot!(
                        dataobject["plot"],
                        [dataobject["previous_x"], dataobject["x"]],
                        [dataobject["previous_y"], dataobject["y"]],
                        # label = "Plot",
                        legend = false,
                        lc = :black,
                        lw = 1,
                    )
                end
                write_pd_mesh(dataobject)
            end
        end
        dataobject["previous_x"] = dataobject["x"]
        dataobject["previous_y"] = dataobject["y"]
        dataobject["previous_extruding"] = true
    end
end

function switch_on(dataobject)
    dataobject["relevant_component"] = true
end
function switch_off(dataobject)
    dataobject["relevant_component"] = false
end
function new_layer(z, dataobject)
    if dataobject["find_min_max"]
        push!(dataobject["layers"], z)
        return
    end
    pd_mesh = dataobject["pd_mesh"]
    if pd_mesh["plot_enabled"]
        dataobject["plot"] = scatter!(
            dataobject["plot"],
            pd_mesh["x_peri"],
            pd_mesh["y_peri"],
            title = "Layer" * string(z),
            xlabel = "X",
            ylabel = "Y",
            ma = 0.5,
            ms = 1,
        )
        savefig(dataobject["plot"], "Output/layer" * string(z) * ".svg")
        dataobject["plot"] = Plots.plot()
        pd_mesh["x_peri"] = []
        pd_mesh["y_peri"] = []
    end

    fill!(pd_mesh["used_ids"], false)
end
function write_pd_mesh(dataobject)
    pd_mesh = dataobject["pd_mesh"]

    pd_mesh["start_point"][1] = dataobject["previous_x"]
    pd_mesh["start_point"][2] = dataobject["previous_y"]
    pd_mesh["point"][1] = dataobject["x"]
    pd_mesh["point"][2] = dataobject["y"]
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])
    v = distance / (dataobject["time"] - dataobject["previous_time"])
    normalize_in_place!(pd_mesh["dir"], pd_mesh["point_diff"])
    neighbors = []
    if distance <= min(pd_mesh["discretization"][1], pd_mesh["discretization"][2])
        idxs = inrange(
            pd_mesh["balltree"],
            pd_mesh["point"],
            max(pd_mesh["discretization"][1], pd_mesh["discretization"][2]),
        )
        # add idxs if distance is less than or equal to discretization
        for i in eachindex(idxs)
            push!(neighbors, idxs[i])
        end
    else
        # @info pd_mesh["start_point"]
        # @info pd_mesh["point"]
        dx =
            pd_mesh["start_point"][1] > pd_mesh["point"][1] ?
            -pd_mesh["discretization"][1] : pd_mesh["discretization"][1]
        dy =
            pd_mesh["start_point"][2] > pd_mesh["point"][2] ?
            -pd_mesh["discretization"][2] : pd_mesh["discretization"][2]
        (xg, yg) = ndgrid(
            pd_mesh["start_point"][1]:dx:pd_mesh["point"][1],
            pd_mesh["start_point"][2]:dy:pd_mesh["point"][2],
        )
        nnodes = length(xg)
        # @info nnodes
        grid = zeros(2, nnodes)
        for i = 1:nnodes
            grid[1, i] = xg[i]
            grid[2, i] = yg[i]
        end
        idxs = inrange(
            pd_mesh["balltree"],
            grid,
            max(pd_mesh["discretization"][1], pd_mesh["discretization"][2]),
        )
        # @info length(idxs)
        # for i in eachindex(xg)
        #     idxs, dists = knn(pd_mesh["balltree"], grid, 4, true)
        for i in eachindex(idxs)
            # @info length(idxs[i])
            for j in eachindex(idxs[i])
                push!(neighbors, idxs[i][j])
            end
        end
        # end
    end
    if length(neighbors) == 0
        return
    end
    # @info neighbors
    for i in eachindex(neighbors)
        if pd_mesh["used_ids"][neighbors[i]]
            continue
        end
        pd_mesh["point"][1] = pd_mesh["grid"][1, neighbors[i]]
        pd_mesh["point"][2] = pd_mesh["grid"][2, neighbors[i]]
        sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
        distance_along_line, distance_to_closest_point = closest_point_to_vector(
            pd_mesh["start_point"],
            pd_mesh["dir"],
            pd_mesh["point"],
            pd_mesh["closest_point"],
            pd_mesh["point_diff"],
        )
        if distance_to_closest_point <= pd_mesh["width"] / 2 &&
           distance_along_line <= distance &&
           distance_along_line >= 0
            time_to_activation = distance_along_line / v
            pd_mesh["used_ids"][neighbors[i]] = true
            push!(
                pd_mesh["mesh_df"],
                [
                    pd_mesh["point"][1],
                    pd_mesh["point"][2],
                    dataobject["z"],
                    1,
                    pd_mesh["volume"],
                    time_to_activation + dataobject["previous_time"],
                ],
            )
            if pd_mesh["plot_enabled"]
                push!(pd_mesh["x_peri"], pd_mesh["point"][1])
                push!(pd_mesh["y_peri"], pd_mesh["point"][2])
            end
        end
    end
end

function main(
    gcode_file::String,
    dx::Float64,
    dy::Float64,
    scale::Float64,
    width::Float64,
    plot_enabled::Bool,
)

    discretization = [dx * scale, dy * scale]

    @info "Read gcode file $gcode_file"
    @info "Params: dx $dx, dy $dy, width $width and scale $scale "
    (grid_x, grid_y), height = write_mesh(gcode_file, true, discretization)

    discretization = [dx * scale, dy * scale, (height / 2) * scale]

    nnodes = length(grid_x)

    grid = zeros(2, nnodes)

    for i = 1:nnodes
        grid[1, i] = grid_x[i]
        grid[2, i] = grid_y[i]
    end

    balltree = BallTree(grid)

    # grid_x = reduce(vcat, xg)
    # grid_y = reduce(vcat, yg)
    # grid_x, grid_y = write_mesh(gcode_file, true, discretization)

    pd_mesh = Dict{String,Any}()
    pd_mesh["plot_enabled"] = plot_enabled
    if plot_enabled
        pd_mesh["x_peri"] = []
        pd_mesh["y_peri"] = []
    end
    pd_mesh["balltree"] = balltree
    pd_mesh["discretization"] = discretization
    pd_mesh["volume"] = discretization[1] * discretization[2] * discretization[3]
    pd_mesh["grid"] = grid
    # pd_mesh["grid_y"] = grid_y
    pd_mesh["previous_time"] = 0
    pd_mesh["previous_x"] = 0
    pd_mesh["previous_y"] = 0
    pd_mesh["previous_z"] = 0
    pd_mesh["previous_extruding"] = 0
    pd_mesh["used_ids"] = fill(false, nnodes)
    pd_mesh["width"] = width

    pd_mesh["mesh_df"] = DataFrame(
        x = Float64[],
        y = Float64[],
        z = Float64[],
        k = Int64[],
        volume = Float64[],
        time = Float64[],
    )
    pd_mesh["closest_point"] = zeros(2)
    pd_mesh["dir"] = zeros(2)
    pd_mesh["start_point"] = zeros(2)
    pd_mesh["point"] = zeros(2)
    pd_mesh["point_diff"] = zeros(2)

    @info "Writing mesh"
    write_mesh(gcode_file, false, discretization, pd_mesh)

    txt_file = joinpath("Output", split(replace(gcode_file, ".gcode" => ".txt"), "/")[end])
    @info "Number of points: $(size(pd_mesh["mesh_df"],1))"
    @info "Printing time: $(maximum(pd_mesh["mesh_df"].time)) seconds"
    write(txt_file, "header: x y z block_id volume Activation_Time\n")
    CSV.write(txt_file, pd_mesh["mesh_df"]; delim = ' ', append = true)

    @info "Finished"
end

parsed_args = parse_commandline()
main(
    parsed_args["filename"],
    parsed_args["dx"],
    parsed_args["dy"],
    parsed_args["scale"],
    parsed_args["width"],
    parsed_args["plot_enabled"],
)
