# SPDX-License-Identifier: BSD-3-Clause
include("./PeriLabInputSchema.jl")
using .Schema
using JSON3

# --- 1. Validation, replacing validate_yaml's inner call ---------------------

example_params = Dict(
    "Blocks" => Dict("block_1" => Dict("Block ID" => 1, "Density" => 2700.0, "Horizon" => 3.0)),
    "Discretization" => Dict("Input Mesh File" => "mesh.txt", "Type" => "Text File"),
    "Solver" => Dict("Verlet" => Dict("Safety Factor" => 0.7),
                     "Final Time" => 1.0, "Initial Time" => 0.0),
)

errors = Schema.validate_params(PERILAB_SCHEMA, example_params)
if isempty(errors)
    println("valid")
else
    foreach(println, errors)
end

# Deliberately broken: two solver kinds at once, missing required Horizon
bad_params = deepcopy(example_params)
bad_params["Solver"]["Static"] = Dict("NLSolve" => true)
delete!(bad_params["Blocks"]["block_1"], "Horizon")

bad_errors = Schema.validate_params(PERILAB_SCHEMA, bad_params)
println("\n--- expected failures ---")
foreach(println, bad_errors)
# Expect: "Solver: must contain only one of: Verlet, Static"
#         "Blocks.block_1.Horizon: missing required key"

# --- 2. JSON Schema export, for PeriHub -------------------------------------

doc = Schema.to_json_schema(PERILAB_SCHEMA, Val(:document);
                            id = "https://perilab.dlr.de/schema/input.json")
open("perilab.schema.json", "w") do io
    JSON3.write(io, doc)
end
println("\nwrote perilab.schema.json")
