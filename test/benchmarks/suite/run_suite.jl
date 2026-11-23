# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

# Benchmark suite: solvers and models of PeriLab on the refined DCB, several MPI rank
# counts and (if available) the distributed matrix backends. Run it once per checkout and
# compare the CSV files with compare.jl:
#
#   julia --project=<checkout> run_suite.jl project=<checkout> label=main \
#         n_points=40000 ranks=1,8,28,56 [cases=verlet_pd_solid,...] [gpu_project=<env>]
#
# Every run is done twice in one process (the first compiles); the timer of the second
# gives t_init (Solver_Manager.init) and t_step (mean of the solver's per-step timer;
# Verlet takes its step count from the stable time step and Final Time). u_norm is
# the norm of the displacements at the end, to check that both checkouts compute the
# same. Cases a checkout does not support (no `Matrix Backend` option) are skipped.

using Exodus
using MPI: mpiexec
using Printf
using YAML
using TOML
include(joinpath(@__DIR__, "create_dcb_mesh.jl"))

args = Dict(split(a, "=") for a in ARGS)
project = get(args, "project", Base.active_project() |> dirname)
label = get(args, "label", basename(project))
n_points = parse(Int64, get(args, "n_points", "40000"))
rank_list = parse.(Int64, split(get(args, "ranks", "1,8,28,56"), ","))
gpu_project = get(args, "gpu_project", "")
gpu_weight = parse(Float64, get(args, "gpu_weight", "55"))
verlet_steps = parse(Int64, get(args, "steps", "50"))
implicit_steps = parse(Int64, get(args, "implicit_steps", "3"))
has_backend = isfile(joinpath(project, "src", "Core", "Solver", "Matrix", "PD_Matrix.jl"))
version = VersionNumber(TOML.parsefile(joinpath(project, "Project.toml"))["version"])

# Per-step timer of a solver as "timer" or "timer|counted by" (v2 times the whole Matrix
# Verlet loop with one timer, the steps are counted by its per-step "solver start").
const V2_TIMERS = Dict("Linear Static Matrix Based" => "Linear Static",
                       "Newmark" => "Newmark step",
                       "Verlet Matrix Based" => "Matrix Verlet|solver start")
step_timer(solver) = version < v"3" ? get(V2_TIMERS, solver, solver) : solver

const ALUMINIUM = Dict("Material Model" => "Correspondence Elastic",
                       "Zero Energy Control" => "Global",
                       "Symmetry" => "isotropic plane stress",
                       "Young's Modulus" => 70000.0, "Poisson's Ratio" => 0.3)
const PD_SOLID = Dict("Material Model" => "PD Solid Elastic",
                      "Symmetry" => "isotropic plane stress",
                      "Young's Modulus" => 70000.0, "Poisson's Ratio" => 0.3)

"""
Case: deck changes on top of the DCB base deck, its rank counts (`nothing` = all) and
whether it needs the distributed matrix backend of this branch.
"""
struct Case
    name::String
    material::Dict
    solver::Dict
    damage::Bool
    thermal::Bool
    serial_only::Bool
    branch_only::Bool
    gpu::String    # "", "GPU", "Mixed"
end
function Case(name, material, solver; damage = false, thermal = false, serial_only = false,
              branch_only = false, gpu = "")
    Case(name, material, solver, damage, thermal, serial_only, branch_only, gpu)
end

verlet() = Dict("Verlet" => Dict("Safety Factor" => 0.9))
# SparseArrays without the new keys, so that main accepts the deck too
function matrix_verlet(backend)
    Dict("Verlet Matrix Based" =>
             backend == "SparseArrays" ?
             Dict("Safety Factor" => 0.9) :
             Dict("Safety Factor" => 0.9,
                  "Matrix Backend" => backend,
                  "GPU Weight" => gpu_weight))
end
function implicit(solver, backend)
    Dict(solver =>
             backend == "SparseArrays" ? Dict() :
             Dict("Matrix Backend" => backend,
                  "Linear Solver" => "BiCGStab",
                  "Linear Solver Tolerance" => 1e-10))
end

const CASES = [Case("verlet_pd_solid", PD_SOLID, verlet()),
    Case("verlet_correspondence", ALUMINIUM, verlet()),
    Case("verlet_pd_solid_damage", PD_SOLID, verlet(); damage = true),
    Case("verlet_pd_solid_thermal", PD_SOLID, verlet(); thermal = true),
    Case("matrix_verlet_sparsearrays", ALUMINIUM, matrix_verlet("SparseArrays");
         serial_only = true),
    Case("linear_static_lu", ALUMINIUM,
         implicit("Linear Static Matrix Based", "SparseArrays");
         serial_only = true),
    Case("newmark_lu", ALUMINIUM, implicit("Newmark", "SparseArrays"); serial_only = true),
    Case("matrix_verlet_cpu", ALUMINIUM, matrix_verlet("CPU"); branch_only = true),
    Case("matrix_verlet_cpu_damage", ALUMINIUM, matrix_verlet("CPU"); damage = true,
         branch_only = true),
    Case("linear_static_bicgstab", ALUMINIUM, implicit("Linear Static Matrix Based", "CPU");
         branch_only = true),
    Case("newmark_bicgstab", ALUMINIUM, implicit("Newmark", "CPU"); branch_only = true),
    Case("matrix_verlet_gpu", ALUMINIUM, matrix_verlet("GPU"); branch_only = true,
         serial_only = true, gpu = "GPU"),
    Case("matrix_verlet_mixed", ALUMINIUM, matrix_verlet("Mixed"); branch_only = true,
         gpu = "Mixed")]

function deck(case::Case, mesh::String, dx::Float64, output::String)
    implicit_solver = !haskey(case.solver, "Verlet") &&
                      !haskey(case.solver, "Verlet Matrix Based")
    block(id) = Dict{String,Any}("Block ID" => id, "Density" => 1.95e-6,
                                 "Horizon" => 3.01dx, "Material Model" => "mat")
    blocks = Dict("block_$i" => block(i) for i in 1:4)
    models = Dict{String,Any}("Material Models" => Dict("mat" => case.material),
                              "Pre Calculation Global" =>
                                  Dict("Bond Associated Deformation Gradient" => false,
                                       "Deformation Gradient" => false,
                                       "Deformed Bond Geometry" => true,
                                       "Shape Tensor" => false))
    if case.damage
        critical = case.material === PD_SOLID ? "Critical Stretch" : "Critical Energy"
        value = critical == "Critical Stretch" ? 1.0e-4 : 0.0035714285714285716
        models["Damage Models"] = Dict("dam" => Dict("Damage Model" => critical,
                                            "Critical Value" => value))
        foreach(i -> blocks["block_$i"]["Damage Model"] = "dam", 1:2)
    end
    if case.thermal
        models["Thermal Models"] = Dict("therm" => Dict("Thermal Model" => "Thermal Flow",
                                             "Type" => "Bond based",
                                             "Thermal Conductivity" => 100.0))
        foreach(b -> (b["Thermal Model"] = "therm"; b["Specific Heat Capacity"] = 900.0),
                values(blocks))
    end
    bcs = Dict("BC_1" => Dict("Coordinate" => "y", "Node Set" => "Node Set 1",
                    "Variable" => "Displacements", "Value" => "4.4*t",
                    "Type" => "Dirichlet"),
               "BC_2" => Dict("Coordinate" => "y", "Node Set" => "Node Set 2",
                    "Variable" => "Displacements", "Value" => "-4.4*t",
                    "Type" => "Dirichlet"))
    case.thermal &&
        (bcs["BC_3"] = Dict("Node Set" => "Node Set 1", "Variable" => "Temperature",
                            "Value" => "100.0", "Type" => "Dirichlet"))
    solver = Dict{String,Any}("Initial Time" => 0.0,
                              "Final Time" => implicit_solver ? 2.0e-4 : 1.0e-5,
                              "Number of Steps" =>
                                  implicit_solver ? implicit_steps :
                                  verlet_steps,
                              "Material Models" => true, "Damage Models" => case.damage,
                              "Thermal Models" => case.thermal)
    merge!(solver, case.solver)
    return Dict("PeriLab" => Dict("Blocks" => blocks, "Boundary Conditions" => bcs,
                     "Models" => models, "Solver" => solver,
                     "Discretization" => Dict("Input Mesh File" => mesh * ".txt",
                          "Type" => "Text File",
                          "Node Sets" => Dict("Node Set 1" => "ns_$(mesh)_1.txt",
                               "Node Set 2" => "ns_$(mesh)_2.txt"),
                          "Bond Filters" =>
                              Dict("bf_1" => Dict("Type" => "Rectangular_Plane",
                                        "Bottom Length" => 2.0,
                                        "Side Length" => 2.0,
                                        "Bottom Unit Vector X" => 1.0,
                                        "Bottom Unit Vector Y" => 0.0,
                                        "Bottom Unit Vector Z" => 0.0,
                                        "Lower Left Corner X" => -0.5,
                                        "Lower Left Corner Y" => 0.0,
                                        "Lower Left Corner Z" => -1.0,
                                        "Normal X" => 0.0,
                                        "Normal Y" => 1.0,
                                        "Normal Z" => 0.0))),
                     "Outputs" =>
                         Dict("Output1" => Dict("Number of Output Steps" => 1,
                                   "Output Filename" => output,
                                   "Output File Type" => "Exodus",
                                   "Output Variables" => Dict("Displacements" => true)))))
end

const WORKER = """
import PeriLab, TimerOutputs, MPI
PeriLab.run(ARGS[1]; silent = true)
PeriLab.run(ARGS[1]; silent = true)
function find(t, name)
    haskey(t.inner_timers, name) && return t.inner_timers[name]
    for c in values(t.inner_timers)
        r = find(c, name)
        r === nothing || return r
    end
    return nothing
end
to = TimerOutputs.get_defaulttimer()
t(name) = TimerOutputs.time(find(to, name)) / 1e9
name, counter = occursin("|", ARGS[3]) ? split(ARGS[3], "|") : (ARGS[3], ARGS[3])
MPI.Comm_rank(MPI.COMM_WORLD) == 0 &&
    write(ARGS[2], string(t("Solver_Manager.init"), " ",
                          t(name) / TimerOutputs.ncalls(find(to, counter))))
"""

function u_norm(file::String)
    exo = ExodusDatabase(file, "r")
    n = length(read_times(exo))
    s = sum(abs2, read_values(exo, NodalVariable, n, 1, "Displacementsx")) +
        sum(abs2, read_values(exo, NodalVariable, n, 1, "Displacementsy"))
    close(exo)
    return sqrt(s)
end

selected = haskey(args, "cases") ? split(args["cases"], ",") : [c.name for c in CASES]
workdir = joinpath(@__DIR__, "runs_$label")
mkpath(workdir)
cd(workdir) do
    mesh, n, dx = create_dcb_mesh(n_points; dir = workdir)
    csv = joinpath(@__DIR__, "results_$(label)_$n.csv")
    isfile(csv) ||
        write(csv, "label,case,ranks,n_points,t_init,t_step,u_norm,status\n")
    for case in CASES
        case.name in selected || continue
        case.branch_only && !has_backend && continue
        !isempty(case.gpu) && isempty(gpu_project) && continue
        env = isempty(case.gpu) ? project : gpu_project
        ranks = case.serial_only ? [1] : rank_list
        case.gpu == "Mixed" && (ranks = filter(>(1), ranks))
        for r in ranks
            output = "$(case.name)_$(r)"
            YAML.write_file(output * ".yaml", deck(case, mesh, dx, output))
            # the per-step timer carries the solver name
            timer = step_timer(only(keys(case.solver)))
            julia = `$(Base.julia_cmd()) --startup-file=no --project=$env -e $WORKER $output.yaml $output.time $timer`
            ok = success(r == 1 ? julia : `$(mpiexec()) -n $r -bind-to core $julia`)
            line = if ok && isfile("$output.time") && isfile("$output.e")
                t_init, t_step = parse.(Float64, split(read("$output.time", String)))
                @sprintf("%s,%s,%d,%d,%.3f,%.6f,%.12e,ok", label, case.name, r, n, t_init,
                         t_step, u_norm(output * ".e"))
            else
                "$label,$(case.name),$r,$n,,,,failed"
            end
            println(line)
            open(io -> println(io, line), csv, "a")
            foreach(rm,
                    filter(f -> startswith(f, output * ".") || startswith(f, output * "_"),
                           readdir()))
        end
    end
end
