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
        "--sampling", "-s"
        help = "sampling"
        arg_type = Float64
        default = 0.2
        "--width", "-w"
        help = "width"
        arg_type = Float64
        default = 0.4
        "--plot_enabled", "-p"
        help = "plot_enabled"
        arg_type = Bool
        default = true
        "--start", "-a"
        help = "start command"
        arg_type = String
        default = nothing
        "--stop", "-o"
        help = "stop command"
        arg_type = String
        default = nothing
        "--end", "-e"
        help = "end command"
        arg_type = String
        default = nothing
        "filename"
        help = "filename"
        required = true
    end

    return parse_args(s)
end

function sub_in_place!(C::Vector{T}, A::Vector{T}, B::Vector{T}) where {T<:Number}
    @assert length(C) == length(A) == length(B)

    @inbounds for i in eachindex(A)
        C[i] = A[i] - B[i]
    end
end

function normalize_in_place!(B::Vector{T}, A::Vector{T}) where {T<:Number}
    nrm = norm(A)
    @inbounds for i in eachindex(A)
        B[i] = A[i] / nrm
    end
end
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
                    pd_mesh = Dict())

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
    callbacks["G4"] = dwell
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
    parseFile(gcode_file, callbacks, myPrinter)

    if pd_mesh["plot_enabled"]
        myPrinter["plot"] = scatter!(myPrinter["plot"],
                                     pd_mesh["x_peri"],
                                     pd_mesh["y_peri"],
                                     title = "Layer" * string(myPrinter["z"]),
                                     xlabel = "X",
                                     ylabel = "Y",
                                     ma = 0.5,
                                     ms = 1)
        savefig(myPrinter["plot"], "Output/layer" * string(myPrinter["z"]) * ".svg")
        myPrinter["plot"] = Plots.plot()
        pd_mesh["x_peri"] = []
        pd_mesh["y_peri"] = []
    end

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
        if dataobject["relevant_component"]
            if dataobject["pd_mesh"]["plot_enabled"]
                dataobject["plot"] = plot!(dataobject["plot"],
                                           [dataobject["previous_x"], dataobject["x"]],
                                           [dataobject["previous_y"], dataobject["y"]],
                                           # label = "Plot",
                                           legend = false,
                                           lc = :black,
                                           lw = 1)
            end
            write_pd_mesh(dataobject)
        end
        dataobject["previous_x"] = dataobject["x"]
        dataobject["previous_y"] = dataobject["y"]
        dataobject["previous_extruding"] = true
    end
end

function dwell(cmds, dataobject)
    s = findfirst((x -> lowercase(x.first) == "s"), cmds)
    p = findfirst((x -> lowercase(x.first) == "p"), cmds)
    wait_time = 0.0
    if s !== nothing
        wait_time = parse(Float64, cmds[s].second)
    end
    if p !== nothing
        wait_time = parse(Float64, cmds[p].second)/1000
    end
    dataobject["previous_time"] = dataobject["time"]
    dataobject["time"] += wait_time
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
    if pd_mesh["plot_enabled"]
        dataobject["plot"] = scatter!(dataobject["plot"],
                                      pd_mesh["x_peri"],
                                      pd_mesh["y_peri"],
                                      title = "Layer" * string(dataobject["z"]),
                                      xlabel = "X",
                                      ylabel = "Y",
                                      ma = 0.5,
                                      ms = 1)
        savefig(dataobject["plot"], "Output/layer" * string(dataobject["z"]) * ".svg")
        dataobject["plot"] = Plots.plot()
        pd_mesh["x_peri"] = []
        pd_mesh["y_peri"] = []
    end
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
    # @info "Rem. Dist1: $(pd_mesh["remaining_distance"])"
    if distance + pd_mesh["remaining_distance"] < pd_mesh["sampling"]
        pd_mesh["remaining_distance"] = pd_mesh["remaining_distance"] + distance
        return
    else
        pd_mesh["remaining_distance"] = pd_mesh["sampling"] - pd_mesh["remaining_distance"]
    end

    # @info "Start1: $(pd_mesh["start_point"])"
    # @info "End: $(pd_mesh["point"])"
    # @info "Distance1: $(distance)"
    # @info "Rem. Dist2: $(pd_mesh["remaining_distance"])"
    # @info "Direction: $(pd_mesh["dir"])"

    pd_mesh["start_point"][1] += pd_mesh["remaining_distance"] * pd_mesh["dir"][1]
    pd_mesh["start_point"][2] += pd_mesh["remaining_distance"] * pd_mesh["dir"][2]
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])

    line_x = []
    line_y = []

    num_of_points_on_line::Int64 = floor(distance / pd_mesh["sampling"]) + 1
    pd_mesh["remaining_distance"] = mod(distance, pd_mesh["sampling"])

    # @info "Start2: $(pd_mesh["start_point"])"
    # @info "Distance2: $(distance)"
    # @info "Point Diff: $(pd_mesh["point_diff"])"
    # @info "Points on line: $(num_of_points_on_line)"
    # @info "Rem. Dist3: $(pd_mesh["remaining_distance"])"

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
        push!(pd_mesh["mesh_df"],
              [
                  pd_mesh["point"][1],
                  pd_mesh["point"][2],
                  dataobject["z"],
                  1,
                  pd_mesh["volume"],
                  time_to_activation + dataobject["previous_time"],
                  180*atan(pd_mesh["dir"][2], pd_mesh["dir"][1])/pi
              ])
        if pd_mesh["plot_enabled"]
            push!(pd_mesh["x_peri"], pd_mesh["point"][1])
            push!(pd_mesh["y_peri"], pd_mesh["point"][2])
        end
    end
end

function main(gcode_file::String,
              sampling::Float64,
              width::Float64,
              plot_enabled::Bool,
              commands_dict)
    @info "Read gcode file $gcode_file"

    pd_mesh = Dict{String,Any}()
    pd_mesh["plot_enabled"] = plot_enabled
    if plot_enabled
        pd_mesh["x_peri"] = []
        pd_mesh["y_peri"] = []
    end
    pd_mesh["sampling"] = sampling
    pd_mesh["volume"] = sampling * width * width
    # pd_mesh["grid"] = grid
    # pd_mesh["grid_y"] = grid_y
    pd_mesh["previous_time"] = 0
    pd_mesh["previous_x"] = 0
    pd_mesh["previous_y"] = 0
    pd_mesh["previous_z"] = 0
    pd_mesh["previous_extruding"] = 0
    pd_mesh["remaining_distance"] = sampling / 2
    pd_mesh["width"] = width

    pd_mesh["mesh_df"] = DataFrame(x = Float64[],
                                   y = Float64[],
                                   z = Float64[],
                                   k = Int64[],
                                   volume = Float64[],
                                   time = Float64[],
                                   angle = Float64[])
    pd_mesh["dir"] = zeros(2)
    pd_mesh["start_point"] = zeros(2)
    pd_mesh["point"] = zeros(2)
    pd_mesh["point_diff"] = zeros(2)

    @info "Writing mesh"
    write_mesh(gcode_file, commands_dict, pd_mesh)

    txt_file = joinpath("Output", split(replace(gcode_file, ".gcode" => ".txt"), "/")[end])
    @info "Number of points: $(size(pd_mesh["mesh_df"],1))"
    @info "Printing time: $(maximum(pd_mesh["mesh_df"].time)) seconds"
    write(txt_file, "header: x y z block_id volume Activation_Time\n")
    CSV.write(txt_file, pd_mesh["mesh_df"]; delim = ' ', append = true)

    @info "Finished"
end

parsed_args = parse_commandline()

commands_dict = Dict{String,Any}()
commands_dict["Start"] = parsed_args["start"]
commands_dict["Stop"] = parsed_args["stop"]
commands_dict["End"] = parsed_args["end"]

main(parsed_args["filename"],
     parsed_args["sampling"],
     parsed_args["width"],
     parsed_args["plot_enabled"], commands_dict)
