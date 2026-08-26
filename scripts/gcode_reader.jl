# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
# Gcode functions taken from MIT Project GcodeParser.jl https://github.com/janvorisek/GcodeParser.jl

using Plots
using LinearAlgebra
using LazyGrids
using Rotations
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
        "--height", "-t"
        help = "height"
        arg_type = Float64
        default = 0.2
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
        "--plot_moves", "-m"
        help = "Plot non-extrusion (travel) moves in a separate color"
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
                callbacks[command](dataObject, z)
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
            # Words that are only ever *parameters* to a motion command, never
            # commands themselves. A line consisting only of these (e.g. "B0
            # C0" right after a G1 move, common in modal g-code that doesn't
            # repeat the G-word on every line) reuses the last motion mode.
            axis_only_letters = ("x", "y", "z", "i", "j", "k", "r", "f", "e",
                                 "a", "b", "c", "u", "v", "w", "s", "p")
            if dataObject !== nothing && lowercase(letter) in axis_only_letters
                command = get(dataObject, "motion_mode", "G0")
            else
                number = cmds[1].second
                if startswith(number, "0")
                    number = string(parse(Int, number))
                end
                command = "$letter$number"
            end
            if haskey(callbacks, command)
                if dataObject === nothing
                    callbacks[command](cmds)
                else
                    callbacks[command](cmds, dataObject)
                    if command in ("G0", "G1", "G2", "G3")
                        dataObject["motion_mode"] = command
                    end
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
    myPrinter["previous_z"] = 0.0
    myPrinter["z"] = 0.0
    myPrinter["b"] = 0.0
    myPrinter["c"] = 0.0
    myPrinter["f"] = 0.0
    myPrinter["e"] = 0.0
    myPrinter["motion_mode"] = "G0"
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
    myPrinter["layer_points"] = Matrix{Int}(undef, 0, 3)
    myPrinter["up_vector"] = [0, 0, 1]

    # Setup a dictionary of callbacks for specified commands
    callbacks = Dict{String,Function}()
    callbacks["G0"] = move # just move the printhead
    callbacks["G1"] = extrude  # move the printhead as well as extrude material
    callbacks["G2"] = arc_cw    # clockwise arc with extrusion
    callbacks["G3"] = arc_ccw   # counter-clockwise arc with extrusion
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

    # Only flush here if there is content that hasn't already been saved by a
    # new_layer() call during parsing (e.g. the file doesn't end on a G0/new
    # layer trigger). Otherwise this would overwrite the last real save with
    # a blank plot, since new_layer() already resets the plot after saving.
    if pd_mesh["plot_enabled"] && !isempty(pd_mesh["x_peri"])
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
    plot_segment!(dataobject, dataobject["previous_x"], dataobject["previous_y"],
                  dataobject["x"], dataobject["y"]; extruding = false)
    dataobject["previous_extruding"] = false
    # dataobject["pd_mesh"]["remaining_distance"] = dataobject["pd_mesh"]["sampling"] / 2
    new_layer(dataobject)
end

# function check_min_max(dataobject, str)
#     if dataobject[str] > dataobject[str * "_max"]
#         dataobject[str * "_max"] = dataobject[str]
#     end
#     if dataobject[str] < dataobject[str * "_min"]
#         dataobject[str * "_min"] = dataobject[str]
#     end
# end

"""
    movement(cmds, dataobject)

Example movement callback for `G0` and `G1` which calculates the total distance moved in all axes.

It is calculated by watching the `X`, `Y` and `Z` axes movement.
"""
function movement(cmds, dataobject; arc_len = nothing)
    dataobject["previous_x"] = dataobject["x"]
    dataobject["previous_y"] = dataobject["y"]
    dataobject["previous_z"] = dataobject["z"]

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
        # if val > dataobject["z"] && x !== nothing
        #     new_layer(val, dataobject)
        # end

        if dataobject["positioning"] === "absolute"
            dz = val - dataobject["z"]
            dataobject["z"] = val
        else
            dz = val
        end
    end

    # B and C are rotary axes (e.g. a tilt/rotation table). This tool has no
    # kinematic model relating them to X/Y/Z, so their values are only
    # tracked for state/reporting purposes and don't contribute to
    # chord_dist or the XY plot.
    b = findfirst((x -> lowercase(x.first) == "b"), cmds)
    if b !== nothing
        val = parse(Float64, cmds[b].second)
        if dataobject["positioning"] === "absolute"
            dataobject["b"] = val
        else
            dataobject["b"] += val
        end
    end

    c = findfirst((x -> lowercase(x.first) == "c"), cmds)
    if c !== nothing
        val = parse(Float64, cmds[c].second)
        if dataobject["positioning"] === "absolute"
            dataobject["c"] = val
        else
            dataobject["c"] += val
        end
    end

    f = findfirst((x -> lowercase(x.first) == "f"), cmds)
    if f !== nothing
        val = parse(Float64, cmds[f].second)
        if val == 0 && x !== nothing
            new_layer(dataobject)
        end

        if dataobject["positioning"] === "absolute"
            dataobject["f"] = val
        end
    end

    chord_dist = sqrt(dx * dx + dy * dy + dz * dz)
    # Use arc_len if provided (for G02/G03 arcs), otherwise chord distance
    move_dist = isnothing(arc_len) ? chord_dist : arc_len
    dataobject["distanceMoved"] += move_dist
    dataobject["previous_time"] = dataobject["time"]
    if dataobject["f"] > 0.0
        dataobject["time"] += move_dist / dataobject["f"] * 60
    end
    return chord_dist
end

"""
    plot_segment!(dataobject, x1, y1, x2, y2; extruding)

Draw one straight XY segment on `dataobject["plot"]`: black when `extruding`
is true, gray when false. No-op if plotting is disabled.
"""
function plot_segment!(dataobject, x1, y1, x2, y2; extruding::Bool)
    pd_mesh = dataobject["pd_mesh"]
    if !pd_mesh["plot_enabled"]
        return
    end
    if extruding
        dataobject["plot"] = plot!(dataobject["plot"], [x1, x2], [y1, y2],
                                   legend = false, lc = :black, lw = 1)
    elseif pd_mesh["plot_moves"]
        dataobject["plot"] = plot!(dataobject["plot"], [x1, x2], [y1, y2],
                                   legend = false, lc = :gray60, lw = 0.8)
    end
end

"""
    extrude(cmds, dataobject)

Example extrusion callback for `G1` which calculates total length of filament extruded.

The extruded filament length is obtained by watching the `E` axis movement in the g-code file.
"""
function extrude(cmds, dataobject)
    distance = movement(cmds, dataobject)

    # calculate used filament length
    e = findfirst((x -> lowercase(x.first) == "e"), cmds)
    is_extruding = false
    de = 0.0
    if e !== nothing
        # Current E axis value
        e_val = parse(Float64, cmds[e].second)

        if dataobject["positioning"] === "absolute"
            de = e_val - dataobject["e"]
            dataobject["e"] = e_val
        else
            de = e_val
        end

        is_extruding = e_val > 0.0
    end

    # A move only actually deposits/removes material if it went somewhere.
    # Pure axis moves (e.g. only B/C rotation, no X/Y/Z change) are not mesh-worthy
    # even if the command is otherwise flagged as "extruding".
    has_motion = distance > 1e-9

    plot_segment!(dataobject, dataobject["previous_x"], dataobject["previous_y"],
                  dataobject["x"], dataobject["y"]; extruding = is_extruding && has_motion)

    if is_extruding && has_motion
        # Used filament
        if e !== nothing
            dataobject["filamentUsage"] += de
        end
        if dataobject["relevant_component"]
            write_pd_mesh(dataobject)
        end
        dataobject["previous_x"] = dataobject["x"]
        dataobject["previous_y"] = dataobject["y"]
        dataobject["previous_z"] = dataobject["z"]
        dataobject["layer_points"] = [dataobject["layer_points"];
                                      [dataobject["x"] dataobject["y"] dataobject["z"]]]
    end
    dataobject["previous_extruding"] = is_extruding && has_motion
end

"""
    arc(cmds, dataobject, clockwise::Bool)

G02/G03 arc interpolation with extrusion. Computes points along a circular arc
defined by the destination (X, Y) and the I, J offset of the arc center.
"""
function arc(cmds, dataobject, clockwise::Bool)
    pd_mesh = dataobject["pd_mesh"]

    # --- Step 1: Parse I/J/R parameters ---
    has_r = false
    r = 0.0
    ic, jc = 0.0, 0.0
    for p in cmds
        lp = lowercase(p.first)
        if lp == "i"
            ic = parse(Float64, p.second)
        elseif lp == "j"
            jc = parse(Float64, p.second)
        elseif lp == "r"
            has_r = true
            r = parse(Float64, p.second)
        end
    end

    # --- Step 2: Pre-compute arc geometry, then call movement() with correct arc length ---
    start_x = dataobject["x"]
    start_y = dataobject["y"]
    start_z = dataobject["z"]

    e_cmd = findfirst((x -> lowercase(x.first) == "e"), cmds)
    has_e = e_cmd !== nothing

    # Parse the destination X/Y directly from the command instead of reading
    # dataobject["x"]/["y"], which still hold the *pre-move* position at this
    # point (movement() hasn't run yet).
    x_cmd = findfirst((p -> lowercase(p.first) == "x"), cmds)
    y_cmd = findfirst((p -> lowercase(p.first) == "y"), cmds)

    if dataobject["positioning"] === "absolute"
        end_x = x_cmd !== nothing ? parse(Float64, cmds[x_cmd].second) : start_x
        end_y = y_cmd !== nothing ? parse(Float64, cmds[y_cmd].second) : start_y
    else
        end_x = start_x + (x_cmd !== nothing ? parse(Float64, cmds[x_cmd].second) : 0.0)
        end_y = start_y + (y_cmd !== nothing ? parse(Float64, cmds[y_cmd].second) : 0.0)
    end

    # Compute center, radius and arc angle before movement — needed to pass correct distance
    if has_r
        # R syntax: compute chord from destination first (R defines the circle; we need end pos for angles)
        dx = end_x - start_x
        dy = end_y - start_y
        chord = sqrt(dx^2 + dy^2)
        if chord < 1e-12
            return
        end
        h = sqrt(max(0.0, r^2 - (chord / 2)^2))
        if clockwise
            cx = (start_x + end_x) / 2 - h * dy / chord
            cy = (start_y + end_y) / 2 + h * dx / chord
        else
            cx = (start_x + end_x) / 2 + h * dy / chord
            cy = (start_y + end_y) / 2 - h * dx / chord
        end
    else
        # I/J center-offset: center known upfront
        cx = start_x + ic
        cy = start_y + jc
        r = sqrt(ic^2 + jc^2)
    end

    θ1 = atan(start_y - cy, start_x - cx)
    θ2 = atan(end_y - cy, end_x - cx)

    # Correct signed angle direction.
    # Clockwise (G02) motion must sweep with *decreasing* angle, so θ2 should
    # end up <= θ1 (subtract 2π if it's currently greater). Counterclockwise
    # (G03) must sweep with *increasing* angle, so θ2 should end up >= θ1
    # (add 2π if it's currently smaller).
    if abs(θ2 - θ1) < 1e-12
        θ2 = θ1 + π * (clockwise ? -1 : 1)
    elseif clockwise
        θ2 > θ1 && (θ2 -= 2π)
    else
        θ2 < θ1 && (θ2 += 2π)
    end

    arc_len = abs(θ2 - θ1) * r

    # Now call movement() with the correct arc length for time/distance tracking
    movement(cmds, dataobject; arc_len = arc_len)

    # --- Step 3: Determine extrusion vs travel ---
    is_extruding = false
    de = 0.0
    if has_e
        e_val = parse(Float64, cmds[e_cmd].second)

        if dataobject["positioning"] === "absolute"
            de = e_val - dataobject["e"]
            dataobject["e"] = e_val
        else
            de = e_val
        end

        is_extruding = e_val > 0.0
    end

    has_motion = arc_len > 1e-9
    is_extruding = is_extruding && has_motion

    if is_extruding && has_e
        dataobject["filamentUsage"] += de
    end

    # --- Step 4: Interpolate along the true curve — every arc is drawn as a
    # curve (gray for travel, black for extrusion), never as a straight
    # secant. Mesh points are only written while extruding AND the current
    # component is relevant (matches extrude()'s behavior).
    if !has_motion
        dataobject["previous_extruding"] = false
        return
    end

    do_mesh = is_extruding && dataobject["relevant_component"]

    n_pts = max(2, floor(Int, arc_len / pd_mesh["sampling"])) + 1
    prev_px = start_x
    prev_py = start_y
    prev_pz = start_z

    for k in 1:(n_pts - 1)
        θ = θ1 + (θ2 - θ1) * (k / (n_pts - 1))
        ax = cx + r * cos(θ)
        ay = cy + r * sin(θ)

        plot_segment!(dataobject, prev_px, prev_py, ax, ay; extruding = is_extruding)
        if do_mesh
            write_pd_mesh_arc(dataobject, prev_px, prev_py, prev_pz, ax, ay,
                              dataobject["z"])
        end
        prev_px = ax
        prev_py = ay
    end

    # Final segment
    plot_segment!(dataobject, prev_px, prev_py, end_x, end_y; extruding = is_extruding)
    if do_mesh
        write_pd_mesh_arc(dataobject, prev_px, prev_py, prev_pz, end_x, end_y,
                          dataobject["z"])

        dataobject["previous_x"] = dataobject["x"]
        dataobject["previous_y"] = dataobject["y"]
        dataobject["previous_z"] = dataobject["z"]
        dataobject["layer_points"] = [dataobject["layer_points"];
                                      [dataobject["x"] dataobject["y"] dataobject["z"]]]
    end
    dataobject["previous_extruding"] = is_extruding
end

arc_cw(cmds, dataobject) = arc(cmds, dataobject, true)
arc_ccw(cmds, dataobject) = arc(cmds, dataobject, false)

function dwell(cmds, dataobject)
    s = findfirst((x -> lowercase(x.first) == "s"), cmds)
    p = findfirst((x -> lowercase(x.first) == "p"), cmds)
    wait_time = 0.0
    if s !== nothing
        wait_time = parse(Float64, cmds[s].second)
    end
    if p !== nothing
        wait_time = parse(Float64, cmds[p].second) / 1000
    end
    dataobject["previous_time"] = dataobject["time"]
    dataobject["time"] += wait_time
    dataobject["previous_extruding"] = false
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
function new_layer(dataobject, z = nothing)
    if dataobject["finsihed"] | !dataobject["relevant_component"]
        return
    end
    pd_mesh = dataobject["pd_mesh"]
    pd_mesh["remaining_distance"] = pd_mesh["sampling"] / 2

    # @info size(dataobject["layer_points"])[1]
    if size(dataobject["layer_points"])[1] != 0
        kdtree = KDTree(transpose(dataobject["layer_points"]))
        point = [dataobject["x"], dataobject["y"], dataobject["z"]]
        if z !== nothing
            point = [dataobject["x"], dataobject["y"], z]
        end
        idx, dist = nn(kdtree, point)
        point_diff = point - dataobject["layer_points"][idx, :]
        dataobject["up_vector"] = point_diff ./ norm(point_diff)
        dataobject["layer_points"] = Matrix{Float64}(undef, 0, 3)
    end

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
"""
    write_pd_mesh_arc(dataobject,
                      sx, sy, sz, ex, ey, ez)

Variant of `write_pd_mesh` for arc segments. Uses the provided start/end
points instead of reading them from `dataobject`.
"""
function write_pd_mesh_arc(dataobject,
                           sx::Number, sy::Number, sz::Number,
                           ex::Number, ey::Number, ez::Number)
    pd_mesh = dataobject["pd_mesh"]

    pd_mesh["start_point"] .= [sx, sy, sz]
    pd_mesh["point"] .= [ex, ey, ez]
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])
    if distance > 1e-12
        roll, pitch,
        yaw = tait_bryant_angles(pd_mesh["point_diff"],
                                 dataobject["up_vector"])
    else
        roll = pitch = yaw = 0.0
    end
    normalize_in_place!(pd_mesh["dir"], pd_mesh["point_diff"])

    if distance + pd_mesh["remaining_distance"] < pd_mesh["sampling"]
        pd_mesh["remaining_distance"] += distance
        return
    else
        pd_mesh["remaining_distance"] = pd_mesh["sampling"] - pd_mesh["remaining_distance"]
    end

    pd_mesh["start_point"][1] += pd_mesh["remaining_distance"] * pd_mesh["dir"][1]
    pd_mesh["start_point"][2] += pd_mesh["remaining_distance"] * pd_mesh["dir"][2]
    pd_mesh["start_point"][3] += pd_mesh["remaining_distance"] * pd_mesh["dir"][3]
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])
    pd_mesh["remaining_distance"] = mod(distance, pd_mesh["sampling"])

    num_of_points::Int64 = floor(distance / pd_mesh["sampling"]) + 1
    if num_of_points > 1
        line_x = collect(range(pd_mesh["start_point"][1],
                               pd_mesh["point"][1] -
                               pd_mesh["remaining_distance"] * pd_mesh["dir"][1],
                               num_of_points))
        line_y = collect(range(pd_mesh["start_point"][2],
                               pd_mesh["point"][2] -
                               pd_mesh["remaining_distance"] * pd_mesh["dir"][2],
                               num_of_points))
        line_z = collect(range(pd_mesh["start_point"][3],
                               pd_mesh["point"][3] -
                               pd_mesh["remaining_distance"] * pd_mesh["dir"][3],
                               num_of_points))
    else
        line_x = [pd_mesh["start_point"][1]]
        line_y = [pd_mesh["start_point"][2]]
        line_z = [pd_mesh["start_point"][3]]
    end

    for i in eachindex(line_x)
        pd_mesh["point"][1] = line_x[i]
        pd_mesh["point"][2] = line_y[i]
        pd_mesh["point"][3] = line_z[i]
        sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
        dist_along_line = distance_along_line(pd_mesh["dir"], pd_mesh["point_diff"])

        v = (dataobject["time"] - dataobject["previous_time"]) > 0 ?
            distance / (dataobject["time"] - dataobject["previous_time"]) : 0.0
        time_to_activation = dist_along_line / v
        push!(pd_mesh["mesh_df"],
              [
                  pd_mesh["point"][1],
                  pd_mesh["point"][2],
                  pd_mesh["point"][3],
                  1,
                  pd_mesh["volume"],
                  time_to_activation + dataobject["previous_time"],
                  roll * 180 / pi,
                  pitch * 180 / pi,
                  yaw * 180 / pi
              ])
        if pd_mesh["plot_enabled"]
            push!(pd_mesh["x_peri"], pd_mesh["point"][1])
            push!(pd_mesh["y_peri"], pd_mesh["point"][2])
        end
    end
end

function tait_bryant_angles(orientation_vector, up_vector = [0, 0, 1])
    forward = orientation_vector ./ norm(orientation_vector)
    right = cross(up_vector, forward)
    right /= norm(right)

    up = cross(forward, right)
    up /= norm(up)

    R = hcat(forward, right, up)
    angles = Rotations.params(RotXYZ(R))

    return angles[1], angles[2], angles[3]
end
function write_pd_mesh(dataobject)
    pd_mesh = dataobject["pd_mesh"]

    pd_mesh["start_point"][1] = dataobject["previous_x"]
    pd_mesh["start_point"][2] = dataobject["previous_y"]
    pd_mesh["start_point"][3] = dataobject["previous_z"]
    pd_mesh["point"][1] = dataobject["x"]
    pd_mesh["point"][2] = dataobject["y"]
    pd_mesh["point"][3] = dataobject["z"]
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])
    roll, pitch, yaw = tait_bryant_angles(pd_mesh["point_diff"], dataobject["up_vector"])
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
    pd_mesh["start_point"][3] += pd_mesh["remaining_distance"] * pd_mesh["dir"][3]
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])

    line_x = []
    line_y = []
    line_z = []

    num_of_points_on_line::Int64 = floor(distance / pd_mesh["sampling"]) + 1
    pd_mesh["remaining_distance"] = mod(distance, pd_mesh["sampling"])

    # @info "Start2: $(pd_mesh["start_point"])"
    # @info "Distance2: $(distance)"
    # @info "Point Diff: $(pd_mesh["point_diff"])"
    # @info "Points on line: $(num_of_points_on_line)"
    # @info "Rem. Dist3: $(pd_mesh["remaining_distance"])"

    if num_of_points_on_line > 1
        line_x = collect(range(pd_mesh["start_point"][1],
                               pd_mesh["point"][1] -
                               pd_mesh["remaining_distance"] * pd_mesh["dir"][1],
                               num_of_points_on_line))
        line_y = collect(range(pd_mesh["start_point"][2],
                               pd_mesh["point"][2] -
                               pd_mesh["remaining_distance"] * pd_mesh["dir"][2],
                               num_of_points_on_line))
        line_z = collect(range(pd_mesh["start_point"][3],
                               pd_mesh["point"][3] -
                               pd_mesh["remaining_distance"] * pd_mesh["dir"][3],
                               num_of_points_on_line))
    else
        line_x = [pd_mesh["start_point"][1]]
        line_y = [pd_mesh["start_point"][2]]
        line_z = [pd_mesh["start_point"][3]]
    end

    if length(line_x) == 0
        return
    end

    for i in eachindex(line_x)
        pd_mesh["point"][1] = line_x[i]
        pd_mesh["point"][2] = line_y[i]
        pd_mesh["point"][3] = line_z[i]
        sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
        dist_along_line = distance_along_line(pd_mesh["dir"], pd_mesh["point_diff"])

        time_to_activation = dist_along_line / v
        push!(pd_mesh["mesh_df"],
              [
                  pd_mesh["point"][1],
                  pd_mesh["point"][2],
                  pd_mesh["point"][3],
                  1,
                  pd_mesh["volume"],
                  time_to_activation + dataobject["previous_time"],
                  roll * 180 / pi,
                  pitch * 180 / pi,
                  yaw * 180 / pi
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
              height::Float64,
              plot_enabled::Bool,
              plot_moves::Bool,
              commands_dict)
    @info "Read gcode file $gcode_file"

    mkpath("Output")

    pd_mesh = Dict{String,Any}()
    pd_mesh["plot_enabled"] = plot_enabled
    pd_mesh["plot_moves"] = plot_moves
    if plot_enabled
        pd_mesh["x_peri"] = []
        pd_mesh["y_peri"] = []
    end
    pd_mesh["sampling"] = sampling
    pd_mesh["volume"] = sampling * width * height
    # pd_mesh["grid"] = grid
    # pd_mesh["grid_y"] = grid_y
    pd_mesh["previous_time"] = 0
    pd_mesh["previous_extruding"] = 0
    pd_mesh["remaining_distance"] = sampling / 2
    pd_mesh["width"] = width
    pd_mesh["height"] = height

    pd_mesh["mesh_df"] = DataFrame(x = Float64[],
                                   y = Float64[],
                                   z = Float64[],
                                   block_id = Int64[],
                                   volume = Float64[],
                                   Activation_Time = Float64[],
                                   Angles_x = Float64[],
                                   Angles_y = Float64[],
                                   Angles_z = Float64[])
    pd_mesh["dir"] = zeros(3)
    pd_mesh["start_point"] = zeros(3)
    pd_mesh["point"] = zeros(3)
    pd_mesh["point_diff"] = zeros(3)

    @info "Writing mesh"
    write_mesh(gcode_file, commands_dict, pd_mesh)

    txt_file = joinpath("Output", split(replace(gcode_file, ".gcode" => ".txt"), "/")[end])
    num_points = size(pd_mesh["mesh_df"], 1)
    @info "Number of points: $(num_points)"
    if num_points == 0
        @info "No points to write. Exiting."
        return
    end
    @info "Printing time: $(maximum(pd_mesh["mesh_df"].Activation_Time)) seconds"
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
     parsed_args["height"],
     parsed_args["plot_enabled"],
     parsed_args["plot_moves"], commands_dict)
