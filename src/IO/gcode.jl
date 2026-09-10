# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
# Gcode functions taken from MIT Project GcodeParser.jl https://github.com/janvorisek/GcodeParser.jl

using LinearAlgebra
using LazyGrids
using Rotations
using CSV, DataFrames
using ProgressBars
using NearestNeighbors
using ..Helpers: sub_in_place!, normalize_in_place!

function distance_along_line(dir::Vector{Float64}, point_diff::Vector{Float64})
    # Calculate the distance from the point to the line segment
    return dot(point_diff, dir) / dot(dir, dir)
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
    iter = progress_bar(0, length(lines) - 1, silent)
    for i in iter
        x = lines[i]
        if occursin(";", x)
            if occursin(";Z:", x)
                command = "new_layer"
                z = parse(Float64, split(x, ":")[2])
                callbacks[command](dataObject, z)
                continue
            else
                command = strip(x[2:end])
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
    callbacks["G1"] = linear  # move the printhead linear
    callbacks["G2"] = arc_cw   # clockwise arc with extrusion
    callbacks["G3"] = arc_ccw  # counter-clockwise arc with extrusion
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
    parseFile(gcode_file, callbacks, myPrinter, silent)

    return
end

# ---------------------------------------------------------------------------
# Shared axis / feedrate / extrusion resolution
# ---------------------------------------------------------------------------

"""
    resolve_axis!(dataobject, cmds, axis::String)

Parse `axis` (one of "x","y","z","b","c") from `cmds` if present and update
`dataobject[axis]` according to the current positioning mode: absolute mode
sets the value directly, relative mode accumulates the delta. Works
uniformly for translational axes (X/Y/Z) and rotary axes (B/C). Returns the
signed delta that was applied (0.0 if the axis word wasn't present).
"""
function resolve_axis!(dataobject, cmds, axis::String)
    idx = findfirst((p -> lowercase(p.first) == axis), cmds)
    idx === nothing && return 0.0

    val = parse(Float64, cmds[idx].second)
    if dataobject["positioning"] === "absolute"
        delta = val - dataobject[axis]
        dataobject[axis] = val
    else
        delta = val
        dataobject[axis] += val
    end
    return delta
end

"""
    resolve_feedrate!(dataobject, cmds)

Parse F from `cmds` if present. Feedrate is always an absolute quantity,
regardless of G90/G91 positioning mode. Preserves the existing convention
that `F0` combined with an `X` word on the same line triggers a new-layer
flush (used by some slicers as a layer-change marker).
"""
function resolve_feedrate!(dataobject, cmds)
    f_idx = findfirst((p -> lowercase(p.first) == "f"), cmds)
    f_idx === nothing && return

    val = parse(Float64, cmds[f_idx].second)
    dataobject["f"] = val

    if val == 0.0 && findfirst((p -> lowercase(p.first) == "x"), cmds) !== nothing
        new_layer(dataobject)
    end
end

"""
    resolve_extrusion!(dataobject, cmds) -> (is_extruding::Bool, de::Float64)

Parse E from `cmds` if present, update `dataobject["e"]`, and report whether
this move deposits material and by how much. `is_extruding` is true whenever
the E word's own value is positive — the relative delta in G91 mode, or the
absolute total in G90 mode. If no `E` word is present, `is_extruding` is
`false` (there is currently no "no-E-in-file means every G1 extrudes"
convention implemented — only explicit positive E values count).
"""
function resolve_extrusion!(dataobject, cmds)
    e_idx = findfirst((p -> lowercase(p.first) == "e"), cmds)
    e_idx === nothing && return false, 0.0

    e_val = parse(Float64, cmds[e_idx].second)
    if dataobject["positioning"] === "absolute"
        de = e_val - dataobject["e"]
        dataobject["e"] = e_val
    else
        de = e_val
        dataobject["e"] += e_val
    end
    return e_val > 0.0, de
end

"""
    advance_clock!(dataobject, path_length)

Updates total distance moved and advances the clock by
`path_length / feedrate`. Must be called *before* any `deposit_mesh_points!`/
`trace_arc!` call for the same move, since mesh deposition reads
`dataobject["time"]`/`["previous_time"]` to compute each point's
Activation_Time — depositing before advancing the clock would timestamp
points using the *previous* move's time window instead of the current one.
"""
function advance_clock!(dataobject, path_length::Float64)
    dataobject["distanceMoved"] += path_length
    dataobject["previous_time"] = dataobject["time"]
    if dataobject["f"] > 0.0
        dataobject["time"] += path_length / dataobject["f"] * 60
    end
end

"""
    finalize_extrusion!(dataobject, is_extruding, path_length, de) -> has_motion::Bool

Called *after* mesh deposition for a move: records whether this move was an
extruding move with actual motion, and — if so — accumulates filament usage
and appends the current position to `layer_points` (used for the next
new_layer's up-vector estimate). Call `advance_clock!` first.
"""
function finalize_extrusion!(dataobject, is_extruding::Bool, path_length::Float64,
                             de::Float64)
    has_motion = path_length > 1e-9
    dataobject["previous_extruding"] = is_extruding && has_motion

    if is_extruding && has_motion
        dataobject["filamentUsage"] += de
        dataobject["layer_points"] = [dataobject["layer_points"];
                                      [dataobject["x"] dataobject["y"] dataobject["z"]]]
    end

    return has_motion
end

# ---------------------------------------------------------------------------
# Mesh deposition (shared by straight moves and each arc sub-chord)
# ---------------------------------------------------------------------------

"""
    deposit_mesh_points!(dataobject, sx, sy, sz, ex, ey, ez)

Sample points along the straight chord from (`sx`,`sy`,`sz`) to
(`ex`,`ey`,`ez`) at `pd_mesh["sampling"]` spacing, carrying over any leftover
distance from the previous call via `pd_mesh["remaining_distance"]`, and push
each sample to `pd_mesh["mesh_df"]` with its estimated activation time,
block id (when `pd_mesh["blocks"]` is set), and orientation angles. Shared by
straight moves and by each small sub-chord of an interpolated arc (see
`trace_arc!`) — replaces the former separate `write_pd_mesh`/
`write_pd_mesh_arc` pair, so block classification and the time-to-activation
guard now behave identically for lines and arcs.
"""
function deposit_mesh_points!(dataobject, sx::Number, sy::Number, sz::Number,
                              ex::Number, ey::Number, ez::Number)
    pd_mesh = dataobject["pd_mesh"]

    pd_mesh["start_point"][1] = sx
    pd_mesh["start_point"][2] = sy
    pd_mesh["start_point"][3] = sz
    pd_mesh["point"][1] = ex
    pd_mesh["point"][2] = ey
    pd_mesh["point"][3] = ez
    sub_in_place!(pd_mesh["point_diff"], pd_mesh["point"], pd_mesh["start_point"])
    distance = norm(pd_mesh["point_diff"])
    if distance < 1e-12
        return
    end

    roll, pitch, yaw = tait_bryant_angles(pd_mesh["point_diff"], dataobject["up_vector"])
    normalize_in_place!(pd_mesh["dir"], pd_mesh["point_diff"])

    dt = dataobject["time"] - dataobject["previous_time"]
    v = dt > 0 ? distance / dt : 0.0

    if distance + pd_mesh["remaining_distance"] < pd_mesh["sampling"]
        pd_mesh["remaining_distance"] += distance
        return
    end
    pd_mesh["remaining_distance"] = pd_mesh["sampling"] - pd_mesh["remaining_distance"]

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
        time_to_activation = v > 0 ? dist_along_line / v : 0.0

        block_id = 1
        if !isnothing(pd_mesh["blocks"])
            global x = pd_mesh["point"][1]
            global y = pd_mesh["point"][2]
            global z = pd_mesh["point"][3]
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
                  pd_mesh["point"][3],
                  block_id,
                  pd_mesh["volume"],
                  time_to_activation + dataobject["previous_time"],
                  roll * 180 / pi,
                  pitch * 180 / pi,
                  yaw * 180 / pi
              ])
    end
end

"""
    trace_arc!(dataobject, cx, cy, r, θ1, θ2,
              start_x, start_y, start_z, end_x, end_y, end_z;
              is_extruding) -> arc_len::Float64

Subdivide the circular arc centered at (`cx`,`cy`) with radius `r`, sweeping
from angle `θ1` to `θ2`, into `pd_mesh["sampling"]`-sized steps, interpolating
Z linearly from `start_z` to `end_z`. Mesh points are only deposited (via
`deposit_mesh_points!`) when `is_extruding` and the current component is
relevant — for a non-extruding arc the subdivision loop is skipped entirely
and the position simply jumps to the endpoint, since there is no plot to
draw. The final substep always lands exactly on (`end_x`,`end_y`,`end_z`) to
avoid trig round-off drift. Returns the total arc length traveled (0.0 for a
degenerate/zero-length arc).
"""
function trace_arc!(dataobject, cx::Float64, cy::Float64, r::Float64,
                    θ1::Float64, θ2::Float64,
                    start_x::Float64, start_y::Float64, start_z::Float64,
                    end_x::Float64, end_y::Float64, end_z::Float64;
                    is_extruding::Bool)
    arc_len = abs(θ2 - θ1) * r
    if arc_len < 1e-9
        dataobject["x"], dataobject["y"], dataobject["z"] = end_x, end_y, end_z
        return 0.0
    end

    pd_mesh = dataobject["pd_mesh"]
    do_mesh = is_extruding && dataobject["relevant_component"]

    if do_mesh
        n_steps = max(1, ceil(Int, arc_len / pd_mesh["sampling"]))
        prev_x, prev_y, prev_z = start_x, start_y, start_z
        for k in 1:n_steps
            if k == n_steps
                nx, ny, nz = end_x, end_y, end_z
            else
                frac = k / n_steps
                θ = θ1 + (θ2 - θ1) * frac
                nx = cx + r * cos(θ)
                ny = cy + r * sin(θ)
                nz = start_z + (end_z - start_z) * frac
            end
            deposit_mesh_points!(dataobject, prev_x, prev_y, prev_z, nx, ny, nz)
            prev_x, prev_y, prev_z = nx, ny, nz
        end
    end

    dataobject["x"], dataobject["y"], dataobject["z"] = end_x, end_y, end_z
    return arc_len
end

# ---------------------------------------------------------------------------
# Motion callbacks
# ---------------------------------------------------------------------------

"""
    move(cmds, dataobject)

G0 callback: non-extruding travel move. Updates X/Y/Z/B/C per the current
positioning mode and triggers a new-layer flush.
"""
function move(cmds, dataobject)
    start_x, start_y, start_z = dataobject["x"], dataobject["y"], dataobject["z"]
    dataobject["previous_x"] = start_x
    dataobject["previous_y"] = start_y
    dataobject["previous_z"] = start_z

    resolve_axis!(dataobject, cmds, "x")
    resolve_axis!(dataobject, cmds, "y")
    resolve_axis!(dataobject, cmds, "z")
    resolve_axis!(dataobject, cmds, "b")
    resolve_axis!(dataobject, cmds, "c")
    resolve_feedrate!(dataobject, cmds)

    path_length = sqrt((dataobject["x"] - start_x)^2 +
                       (dataobject["y"] - start_y)^2 +
                       (dataobject["z"] - start_z)^2)

    advance_clock!(dataobject, path_length)
    finalize_extrusion!(dataobject, false, path_length, 0.0)
    new_layer(dataobject)
end

"""
    linear(cmds, dataobject)

G1 callback. Two cases:
- If the line contains a `C` word, it's treated as a cylindrical rotation
  about the origin: the current radius (distance from origin) is held fixed
  while C sweeps from its old angle to its new one, traced via `trace_arc!`
  so the deposited path follows the true arc rather than a straight chord.
- Otherwise it's a straight-line move between the old and new X/Y/Z.

Extrusion is governed by an `E` word (see `resolve_extrusion!`). The clock is
advanced via `advance_clock!` *before* any mesh deposition (so activation
times are computed against this move's own time window), and filament/
layer_points bookkeeping is finalized afterward via `finalize_extrusion!`.
"""
function linear(cmds, dataobject)
    start_x, start_y, start_z = dataobject["x"], dataobject["y"], dataobject["z"]
    dataobject["previous_x"] = start_x
    dataobject["previous_y"] = start_y
    dataobject["previous_z"] = start_z

    has_c = findfirst((p -> lowercase(p.first) == "c"), cmds) !== nothing
    old_c = dataobject["c"]

    resolve_axis!(dataobject, cmds, "x")
    resolve_axis!(dataobject, cmds, "y")
    resolve_axis!(dataobject, cmds, "z")
    resolve_axis!(dataobject, cmds, "b")
    resolve_axis!(dataobject, cmds, "c")
    resolve_feedrate!(dataobject, cmds)

    is_extruding, de = resolve_extrusion!(dataobject, cmds)

    # First determine the geometry (endpoint + path length) without touching
    # the clock or the mesh yet. (if/else doesn't introduce its own scope in
    # Julia, so r/θ1/θ2/end_x/end_y/end_z/path_length assigned below remain
    # visible for the rest of the function.)
    if has_c
        # Cylindrical mapping: C rotates the current radius about the
        # origin. Radius and start Z are taken from the position *before*
        # this line's updates; the new C value (already resolved above)
        # gives the sweep's end angle.
        r = sqrt(start_x^2 + start_y^2)
        θ1 = deg2rad(old_c)
        θ2 = deg2rad(dataobject["c"])
        end_x = r * cos(θ2)
        end_y = r * sin(θ2)
        end_z = dataobject["z"]
        path_length = abs(θ2 - θ1) * r
    else
        end_x, end_y, end_z = dataobject["x"], dataobject["y"], dataobject["z"]
        dx = end_x - start_x
        dy = end_y - start_y
        dz = end_z - start_z
        path_length = sqrt(dx^2 + dy^2 + dz^2)
    end

    # Advance the clock BEFORE depositing mesh points: deposit_mesh_points!
    # reads dataobject["time"]/["previous_time"] to compute each point's
    # Activation_Time, so the clock must already reflect *this* move.
    advance_clock!(dataobject, path_length)

    if has_c
        trace_arc!(dataobject, 0.0, 0.0, r, θ1, θ2,
                   start_x, start_y, start_z,
                   end_x, end_y, end_z; is_extruding = is_extruding)
    elseif is_extruding && path_length > 1e-9 && dataobject["relevant_component"]
        deposit_mesh_points!(dataobject, start_x, start_y, start_z,
                             end_x, end_y, end_z)
    end

    finalize_extrusion!(dataobject, is_extruding, path_length, de)
end

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
end

"""
    arc(cmds, dataobject, clockwise::Bool)

G02/G03 arc interpolation with extrusion. Computes the arc center and radius
from either I/J (center offset from the start point) or R (radius, with
center chosen to match the requested rotation direction), then traces it via
`trace_arc!`. Extrusion is governed by an `E` word exactly as in `linear`.
"""
function arc(cmds, dataobject, clockwise::Bool)
    start_x, start_y, start_z = dataobject["x"], dataobject["y"], dataobject["z"]
    dataobject["previous_x"] = start_x
    dataobject["previous_y"] = start_y
    dataobject["previous_z"] = start_z

    # --- Step 1: Parse I/J/R parameters ---
    has_r = false
    r_param = 0.0
    ic, jc = 0.0, 0.0
    for p in cmds
        lp = lowercase(p.first)
        if lp == "i"
            ic = parse(Float64, p.second)
        elseif lp == "j"
            jc = parse(Float64, p.second)
        elseif lp == "r"
            has_r = true
            r_param = parse(Float64, p.second)
        end
    end

    resolve_axis!(dataobject, cmds, "x")
    resolve_axis!(dataobject, cmds, "y")
    resolve_axis!(dataobject, cmds, "z")
    resolve_axis!(dataobject, cmds, "b")
    resolve_axis!(dataobject, cmds, "c")
    resolve_feedrate!(dataobject, cmds)

    end_x, end_y, end_z = dataobject["x"], dataobject["y"], dataobject["z"]

    # --- Step 2: Compute center and radius ---
    if has_r
        dx = end_x - start_x
        dy = end_y - start_y
        chord = sqrt(dx^2 + dy^2)
        if chord < 1e-12
            dataobject["previous_extruding"] = false
            return  # zero-length chord, cannot compute arc
        end
        h = sqrt(max(0.0, r_param^2 - (chord / 2)^2))
        if clockwise
            cx = (start_x + end_x) / 2 - h * dy / chord
            cy = (start_y + end_y) / 2 + h * dx / chord
        else
            cx = (start_x + end_x) / 2 + h * dy / chord
            cy = (start_y + end_y) / 2 - h * dx / chord
        end
        r = r_param
    else
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
        θ2 = θ1 + π * (clockwise ? -1.0 : 1.0)
    elseif clockwise
        θ2 > θ1 && (θ2 -= 2π)
    else
        θ2 < θ1 && (θ2 += 2π)
    end

    is_extruding, de = resolve_extrusion!(dataobject, cmds)

    arc_len = abs(θ2 - θ1) * r

    # Advance the clock BEFORE trace_arc! deposits mesh points (see note in
    # linear()) — otherwise every arc point's Activation_Time is computed
    # against the previous move's time window instead of this one's.
    advance_clock!(dataobject, arc_len)

    path_length = trace_arc!(dataobject, cx, cy, r, θ1, θ2,
                             start_x, start_y, start_z,
                             end_x, end_y, end_z; is_extruding = is_extruding)

    finalize_extrusion!(dataobject, is_extruding, path_length, de)
end

arc_cw(cmds, dataobject) = arc(cmds, dataobject, true)
arc_ccw(cmds, dataobject) = arc(cmds, dataobject, false)

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
function new_layer(dataobject, z = -1)
    if dataobject["finsihed"] | !dataobject["relevant_component"]
        return
    end
    pd_mesh = dataobject["pd_mesh"]
    pd_mesh["remaining_distance"] = pd_mesh["sampling"] / 2

    if size(dataobject["layer_points"])[1] != 0
        kdtree = KDTree(transpose(dataobject["layer_points"]))
        point = [dataobject["x"], dataobject["y"], dataobject["z"]]
        if z != -1
            point = [dataobject["x"], dataobject["y"], z]
        end
        idx, dist = nn(kdtree, point)
        point_diff = point - dataobject["layer_points"][idx, :]
        dataobject["up_vector"] = point_diff ./ norm(point_diff)
        dataobject["layer_points"] = Matrix{Int}(undef, 0, 3)
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

function get_gcode_mesh(gcode_file::String, params::Dict, silent)
    sampling = params["Discretization"]["Gcode"]["Sampling"]
    scale = get(params["Discretization"]["Gcode"], "Scale", 1)
    width = params["Discretization"]["Gcode"]["Width"]
    height = params["Discretization"]["Gcode"]["Height"]
    blocks = get(params["Discretization"]["Gcode"], "Blocks", nothing)

    commands_dict = Dict{String,Any}()
    commands_dict["Start"] = get(params["Discretization"]["Gcode"], "Start Command",
                                 nothing)
    commands_dict["Stop"] = get(params["Discretization"]["Gcode"], "Stop Command", nothing)
    commands_dict["End"] = get(params["Discretization"]["Gcode"], "End Command", nothing)

    if !isnothing(commands_dict["Start"])
        if isnothing(commands_dict["Stop"])
            @abort "Start command is set but no stop command"
        end
    end
    if !isnothing(commands_dict["Stop"])
        if isnothing(commands_dict["Start"])
            @abort "Stop command is set but no start command"
        end
    end

    @info "Read gcode file $gcode_file"
    @info "Params: Sampling $sampling, width $width and scale $scale "

    pd_mesh = Dict{String,Any}()
    pd_mesh["sampling"] = sampling
    pd_mesh["volume"] = sampling * width * height
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
    pd_mesh["dir"] = zeros(3)
    pd_mesh["start_point"] = zeros(3)
    pd_mesh["point"] = zeros(3)
    pd_mesh["point_diff"] = zeros(3)

    @info "Writing mesh"
    write_mesh(gcode_file, commands_dict, silent, pd_mesh)

    if size(pd_mesh["mesh_df"], 1) == 0
        @abort "No points found in the gcode file, maybe the gcode format is not supported?"
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
