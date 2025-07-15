# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
# Gcode functions taken from MIT Project GcodeParser.jl https://github.com/janvorisek/GcodeParser.jl

using LinearAlgebra
using LazyGrids
using CSV, DataFrames
using ProgressBars
using NearestNeighbors
include("../Support/Helpers.jl")
using .Helpers: sub_in_place!, normalize_in_place!

function distance_along_line(dir::Vector{Float64}, point_diff::Vector{Float64})

    # Calculate the distance from the point to the line segment
    distance_along_line = dot(point_diff, dir) / dot(dir, dir)

    return distance_along_line
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
function parseLine(line::String,
                   returnPair::Bool = true)::Array{Union{String,Pair{String,String}},1}
    line = stripComments(line)

    # Match commands
    gcode_regex = r"/(%.*)|({.*)|((?:\$\$)|(?:\$[a-zA-Z0-9#]*))|([a-zA-Z][0-9\+\-\.]+)|(\*[0-9]+)/igm"

    # array of matched strings
    matches = collect(String(m.match) for m in eachmatch(gcode_regex, line))

    if returnPair
        return collect(first(m, 1) => last(m, length(m) - 1) for m in matches)
    end

    return matches
end

function parseFile(path::String, callbacks::Dict{String,Function}, dataObject, silent)
    lines = readlines(path)
    # open(path) do f
    #     line = 1
    #     while !eof(f)
    #         x = readline(f);
    iter = progress_bar(0, length(lines)-1, silent)
    for i in iter
        x = lines[i]
        if occursin(";", x)
            if occursin(";Z:", x)
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

function write_mesh(gcode_file, commands_dict,
                    silent = false, pd_mesh = Dict())

    # create any data object
    # it will be passed as a second parameter to your callbacks
    # here simple dictionary is used to store information during the print
    myPrinter = Dict{String,Any}()
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
    myPrinter["relevant_component"] = true
    myPrinter["x_min"] = 1.e100
    myPrinter["x_max"] = 0.0
    myPrinter["y_min"] = 1.e100
    myPrinter["y_max"] = 0.0
    myPrinter["pd_mesh"] = pd_mesh
    myPrinter["layers"] = []
    myPrinter["finsihed"] = false

    # Setup a dictionary of callbacks for specified commands
    callbacks = Dict{String,Function}()
    callbacks["G0"] = move # just move the printhead
    callbacks["G1"] = extrude  # move the printhead as well as extrude material
    callbacks["new_layer"] = new_layer

    if !isnothing(commands_dict["Start"])
        for command in split(commands_dict["Start"], ",")
            callbacks[command] = switch_on
        end
        myPrinter["relevant_component"] = false
    end
    if !isnothing(commands_dict["Stop"])
        for command in split(commands_dict["Stop"], ",")
            callbacks[command] = switch_off
        end
    end
    if !isnothing(commands_dict["End"])
        callbacks[commands_dict["End"]] = finished
    end

    # watch out for relative and absolute positioning
    callbacks["G90"] = (cmds, dataobject) -> dataobject["positioning"] = "absolute"
    callbacks["G91"] = (cmds, dataobject) -> dataobject["positioning"] = "relative"

    # parse g-code file and simulate print using our own callbacks and data object
    parseFile(gcode_file, callbacks, myPrinter, silent)

    return
end

function move(cmds, dataobject)
    movement(cmds, dataobject)
    dataobject["previous_extruding"] = false
    dataobject["pd_mesh"]["remaining_distance"] = dataobject["pd_mesh"]["sampling"] / 2
end

function check_min_max(dataobject, str)
    if dataobject[str] > dataobject[str * "_max"]
        dataobject[str * "_max"] = dataobject[str]
    end
    if dataobject[str] < dataobject[str * "_min"]
        dataobject[str * "_min"] = dataobject[str]
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
        if val > dataobject["z"] && x !== nothing
            new_layer(val, dataobject)
        end

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
    if dataobject["f"] > 0.0
        dataobject["time"] += distance / dataobject["f"] * 60
    end
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
        if dataobject["relevant_component"] # && dataobject["previous_extruding"]
            write_pd_mesh(dataobject)
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
function finished(dataobject)
    dataobject["relevant_component"] = false
    dataobject["finsihed"] = true
end
function new_layer(z, dataobject)
    if dataobject["finsihed"]
        return
    end
    pd_mesh = dataobject["pd_mesh"]
    pd_mesh["remaining_distance"] = pd_mesh["sampling"] / 2
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
    if distance + pd_mesh["remaining_distance"] < pd_mesh["sampling"]
        pd_mesh["remaining_distance"] = pd_mesh["remaining_distance"] + distance
        return
    else
        pd_mesh["remaining_distance"] = pd_mesh["sampling"] - pd_mesh["remaining_distance"]
    end

    pd_mesh["start_point"][1] += pd_mesh["remaining_distance"] * pd_mesh["dir"][1]
    pd_mesh["start_point"][2] += pd_mesh["remaining_distance"] * pd_mesh["dir"][2]
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])

    line_x = []
    line_y = []

    num_of_points_on_line::Int64 = floor(distance / pd_mesh["sampling"]) + 1
    pd_mesh["remaining_distance"] = mod(distance, pd_mesh["sampling"])

    if num_of_points_on_line > 1
        line_x = collect(range(pd_mesh["start_point"][1],
                               pd_mesh["point"][1]-pd_mesh["remaining_distance"]*pd_mesh["dir"][1],
                               num_of_points_on_line))
        line_y = collect(range(pd_mesh["start_point"][2],
                               pd_mesh["point"][2]-pd_mesh["remaining_distance"]*pd_mesh["dir"][2],
                               num_of_points_on_line))
    else
        line_x = [pd_mesh["start_point"][1]]
        line_y = [pd_mesh["start_point"][2]]
    end

    if length(line_x) == 0
        return
    end

    for i in eachindex(line_x)
        pd_mesh["point"][1] = line_x[i]
        pd_mesh["point"][2] = line_y[i]
        sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
        dist_along_line = distance_along_line(pd_mesh["dir"], pd_mesh["point_diff"])

        time_to_activation = dist_along_line / v
        block_id = 1
        if !isnothing(pd_mesh["blocks"])
            global x = pd_mesh["point"][1]
            global y = pd_mesh["point"][2]
            global z = dataobject["z"]
            for block in pd_mesh["blocks"]
                if eval(Meta.parse(block[2]))
                    block_id = block[1]
                end
            end
        end
        push!(pd_mesh["mesh_df"],
              [
                  pd_mesh["point"][1],
                  pd_mesh["point"][2],
                  dataobject["z"],
                  block_id,
                  pd_mesh["volume"],
                  time_to_activation + dataobject["previous_time"],
                  0.0,
                  0.0,
                  180*atan(pd_mesh["dir"][2], pd_mesh["dir"][1])/pi
              ])
    end
end

function get_gcode_mesh(gcode_file::String, params::Dict, silent)
    sampling = params["Discretization"]["Gcode"]["Sampling"]
    scale = get(params["Discretization"]["Gcode"], "Scale", 1)
    width = params["Discretization"]["Gcode"]["Width"]
    blocks = get(params["Discretization"]["Gcode"], "Blocks", nothing)

    commands_dict = Dict{String,Any}()
    commands_dict["Start"] = get(params["Discretization"]["Gcode"], "Start Command",
                                 nothing)
    commands_dict["Stop"] = get(params["Discretization"]["Gcode"], "Stop Command", nothing)
    commands_dict["End"] = get(params["Discretization"]["Gcode"], "End Command", nothing)

    if !isnothing(commands_dict["Start"])
        if isnothing(commands_dict["Stop"])
            @error "Start command is set but no stop command"
        end
    end
    if !isnothing(commands_dict["Stop"])
        if isnothing(commands_dict["Start"])
            @error "Stop command is set but no start command"
        end
    end

    @info "Read gcode file $gcode_file"
    @info "Params: Sampling $sampling, width $width and scale $scale "

    pd_mesh = Dict{String,Any}()
    pd_mesh["sampling"] = sampling
    pd_mesh["volume"] = sampling * width * width
    pd_mesh["previous_time"] = 0
    pd_mesh["previous_x"] = 0
    pd_mesh["previous_y"] = 0
    pd_mesh["previous_z"] = 0
    pd_mesh["previous_extruding"] = 0
    pd_mesh["width"] = width
    pd_mesh["remaining_distance"] = sampling / 2
    pd_mesh["blocks"] = blocks

    pd_mesh["mesh_df"] = DataFrame(x = Float64[],
                                   y = Float64[],
                                   z = Float64[],
                                   block_id = Int64[],
                                   volume = Float64[],
                                   Activation_Time = Float64[],
                                   Angles_x = Float64[],
                                   Angles_y = Float64[],
                                   Angles_z = Float64[])
    pd_mesh["dir"] = zeros(2)
    pd_mesh["start_point"] = zeros(2)
    pd_mesh["point"] = zeros(2)
    pd_mesh["point_diff"] = zeros(2)

    @info "Writing mesh"
    write_mesh(gcode_file, commands_dict, silent, pd_mesh)

    if size(pd_mesh["mesh_df"], 1) == 0
        @error "No points found in the gcode file, maybe the gcode format is not supported?"
        return nothing
    end
    @info "Number of points: $(size(pd_mesh["mesh_df"],1))"
    @info "Printing time: $(maximum(pd_mesh["mesh_df"].Activation_Time)) seconds"

    txt_file = replace(gcode_file, ".gcode" => ".txt")
    write(txt_file,
          "header: x y z block_id volume Activation_Time Angles_x Angles_y Angles_z\n")
    CSV.write(txt_file, pd_mesh["mesh_df"]; delim = ' ', append = true)

    @info "Finished reading mesh data"

    return pd_mesh["mesh_df"]
end
