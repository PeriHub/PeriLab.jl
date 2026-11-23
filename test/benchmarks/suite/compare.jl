# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

# Markdown table of run_suite.jl results, the first file is the baseline of the speedups:
#   julia compare.jl results_v2.2.5_N.csv results_main_N.csv results_branch_N.csv
# u_norm: "=" if all checkouts that ran the case agree to 1e-6, else the values.
using Printf

runs = [[split(l, ",") for l in readlines(f)[2:end]] for f in ARGS]
labels = [r[1][1] for r in runs]
key(r) = (r[2], parse(Int, r[3]))
index = [Dict(key(r) => r for r in rs) for rs in runs]
keys_ordered = unique(key(r) for rs in reverse(runs) for r in rs)
num(x) = isempty(x) ? NaN : parse(Float64, x)
cell(x, fmt) = isnan(x) ? "–" : Printf.format(Printf.Format(fmt), x)

println("| case | ranks | ", join(["t_step $l [s]" for l in labels], " | "), " | ",
        join(["speedup $l" for l in labels[2:end]], " | "), " | u_norm |")
println("|---|---:|", repeat("---:|", 2 * length(labels) - 1), "---|")
for k in keys_ordered
    rows = [get(ix, k, nothing) for ix in index]
    t = [r === nothing ? NaN : num(r[6]) for r in rows]
    u = [r === nothing ? NaN : num(r[7]) for r in rows]
    present = filter(!isnan, u)
    same = all(x -> isapprox(x, present[end]; rtol = 1e-6), present)
    u_cell = same ? "=" : join([cell(x, "%.3e") for x in u], " / ")
    println("| $(k[1]) | $(k[2]) | ", join([cell(x, "%.5f") for x in t], " | "), " | ",
            join([cell(t[1] / x, "%.2f") for x in t[2:end]], " | "), " | $u_cell |")
end
