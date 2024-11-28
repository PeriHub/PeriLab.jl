# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("standard.jl")
using Plots
mesh = get_mesh_header()
node = get_node_header()

x = 0.2 # m
y = 0.2 # m
notch_x = 0.025
notch_y = 0.005
dx = 0.002
vol = dx * dx


# Parameter
x_min, x_max = -0.1, 0.1  # Grenzen in x-Richtung
y_min, y_max = -0.1, 0.1  # Grenzen in y-Richtung


# Gitterpunkte berechnen
x_vals = x_min:dx:x_max
y_vals = y_min:dx:y_max

# Alle Kombinationen der (x, y)-Werte in eine Liste speichern
mesh_entry = [(x, y, 1, vol) for x in x_vals for y in y_vals]


# Notches
notch_left = (-0.1, -0.0025, -0.075, 0.0025)  # (x_min, y_min, x_max, y_max) links
notch_right = (0.075, -0.0025, 0.1, 0.0025)  # (x_min, y_min, x_max, y_max) rechts
function in_notch(x, y, notch)
    x_min, y_min, x_max, y_max = notch
    return x_min ≤ x ≤ x_max && y_min ≤ y ≤ y_max
end

mesh_entry = [
    (x, y, 1, vol) for x in x_vals for
    y in y_vals if !in_notch(x, y, notch_left) && !in_notch(x, y, notch_right)
]
#scatter([p[1] for p in punkte], [p[2] for p in mesh_entry], legend=false, xlabel="x", ylabel="y")


open("mesh_notch.txt", "w") do io
    write(io, mesh)  # Header
    for (x, y, blockID, volume) in mesh_entry
        write(io, "$x $y, $blockID, $volume\n")
    end
end
function point_ids(points, condition)
    return [i for (i, (x, y)) in enumerate(points) if condition(x, y)]
end
# 1. Punkte auf einer Linie (y = wert)
y_val = -0.1
ids_1 = point_ids(mesh_entry, (x, y) -> y == y_val)

y_val = 0.1
ids_2 = point_ids(mesh_entry, (x, y) -> y == y_val)

x_const = -0.1
y_val = 0
ids_3 = point_ids(mesh_entry, (x, y) -> x == x_const && y > y_val)

x_const = 0.1
y_val = 0
ids_4 = point_ids(mesh_entry, (x, y) -> x == x_const && y < y_val)

open("nodes_1.txt", "w") do io
    write(io, node)  # Header
    for id in ids_1
        write(io, "$id\n")
    end
end

open("nodes_2.txt", "w") do io
    write(io, node)  # Header
    for id in ids_2
        write(io, "$id\n")
    end
end
open("nodes_3.txt", "w") do io
    write(io, node)  # Header
    for id in ids_3
        write(io, "$id\n")
    end
end
open("nodes_4.txt", "w") do io
    write(io, node)  # Header
    for id in ids_4
        write(io, "$id\n")
    end
end
