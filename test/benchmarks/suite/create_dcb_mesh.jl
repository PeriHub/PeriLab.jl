# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

# Refined DCB of test_linear_static_matrix_based_solver/DCBmodel.txt (same outer geometry,
# blocks, node sets and bond filter, grid spacing scaled to about `n_points` points).
# julia create_dcb_mesh.jl n_points=100000  →  DCB_<n>.txt, ns_DCB_<n>_1.txt, ns_DCB_<n>_2.txt

using Printf

"""
    create_dcb_mesh(n_points::Int64; dir = @__DIR__) -> (name, n, dx)

Grid of `21 s × 24 s` points (`s = √(n_points / 504)`): x ∈ [0, 21 dx₀), two halves
`0 < |y| < 12 dx₀` with the gap of the original model at y = 0, `dx₀ = 3/11`. Blocks 3
(top) and 4 (bottom) are the loading arms `x < 3 dx₀`, node set 1 = block 3, node set 2 =
block 4. Volume `dx²`; the horizon of the deck is `3.01 dx`.
"""
function create_dcb_mesh(n_points::Int64; dir::String = @__DIR__)
    dx0 = 3 / 11
    s = sqrt(n_points / 504)
    nx = round(Int64, 21s)
    ny = round(Int64, 12s)          # rows per half
    dx = 21dx0 / nx
    name = "DCB_$(2nx * ny)"
    ns = (Int64[], Int64[])
    open(joinpath(dir, name * ".txt"), "w") do io
        println(io, "header: x y block_id volume")
        id = 0
        for half in (-1, 1), j in 1:ny, i in 1:nx
            x = (i - 1) * dx
            y = half * (j - 0.5) * dx
            arm = x < 3dx0 - 1e-12
            block = arm ? (half > 0 ? 3 : 4) : (half > 0 ? 2 : 1)
            id += 1
            arm && push!(ns[half > 0 ? 1 : 2], id)
            @printf(io, "%.15e %.15e %d %.15e\n", x, y, block, dx^2)
        end
    end
    for k in 1:2
        open(joinpath(dir, "ns_$(name)_$k.txt"), "w") do io
            println(io, "header: global_id")
            foreach(i -> println(io, i), ns[k])
        end
    end
    return name, 2nx * ny, dx
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = Dict(split(a, "=") for a in ARGS)
    name, n, dx = create_dcb_mesh(parse(Int64, get(args, "n_points", "100000")))
    println("$name: $n points, dx = $dx")
end
