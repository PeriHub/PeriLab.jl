# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    PeriLab

A module for managing and executing peridynamic simulations in PeriLab.

This module provides functionality for running simulations in the PeriLab environment. It includes functions for initializing simulations, processing data, running solvers, and managing results.

## Modules

- `Core/data_manager.jl`: Data manager module for data management and access.
- `IO/logging.jl`: Module for setting up and managing logging.
- `IO/IO.jl`: Input/output functions for handling data files.
- `Core/Solver/Solver_control.jl`: Solver control module for managing simulation solvers.

## Dependencies

This module depends on the following external packages and modules:
- `MPI`: Message Passing Interface for distributed computing.
- `Pkg`: Julia's package manager for managing project dependencies.
- `TimerOutputs`: Module for measuring and displaying code execution times.
- `Logging`: Julia's built-in logging framework.
- `ArgParse`: Module for parsing command-line arguments.

## Usage

To run a simulation using PeriLab, you can use the `main()` function, which takes several optional parameters to control the simulation process. You can specify whether the simulation should be a dry run, enable verbose output, enable debugging, and run in silent mode.

For example:
```julia
main("examples/Dogbone/Dogbone.yaml"; output_dir="", dry_run=false, verbose=false, debug=false, silent=false, reload=false)
"""

module PeriLab
include("./Core/data_manager.jl")
include("./IO/logging.jl")
include("./IO/IO.jl")
include("./Core/Solver/Solver_control.jl")

using MPI
using TimerOutputs
using Logging
using ArgParse
using Dates
using LibGit2

const to = TimerOutput()
using .Data_manager

# import PrecompileTools
import .Logging_module
import .IO
import .Solver

PERILAB_VERSION = "1.2.2"

export main

"""
    print_banner()

Prints a banner displaying information about the PeriLab application.

This function prints a banner containing details about the PeriLab application, including its name, version, copyright, contact information, and license. It provides a visual introduction to the application.
"""
function print_banner()
    println(
        """\e[]
\e[1;36mPeriLab. \e[0m                  \e[1;32md8b \e[1;36m888               888\e[0m       |  Version: $PERILAB_VERSION
\e[1;36m888   Y88b\e[0m                 \e[1;32mY8P \e[1;36m888               888\e[0m       |  Copyright:
\e[1;36m888    888\e[0m                     \e[1;36m888               888\e[0m       |  Dr.-Ing. Christian Willberg (https://orcid.org/0000-0003-2433-9183)
\e[1;36m888   d88P\e[0m \e[1;36m.d88b.\e[0m  \e[1;36m888d888 888 888       \e[1;36m8888b.\e[0m  \e[1;36m88888b.\e[0m   |  M.Sc. Jan-Timo Hesse (https://orcid.org/0000-0002-3006-1520)
\e[1;36m8888888P"\e[0m \e[1;36md8P  Y8b\e[0m \e[1;36m888P"   888 888          \e[1;36m"88b\e[0m \e[1;36m888 "88b\e[0m  |  Contact: christian.willberg@h2.de, jan-timo.hesse@dlr.de
\e[1;36m888\e[0m       \e[1;36m88888888\e[0m \e[1;36m888\e[0m     \e[1;36m888\e[0m \e[1;36m888\e[0m      \e[1;36m.d888888\e[0m \e[1;36m888  888\e[0m  |  License: BSD-3-Clause
\e[1;36m888\e[0m       \e[1;36mY8b.\e[0m     \e[1;36m888\e[0m     \e[1;36m888\e[0m \e[1;36m888\e[0m      \e[1;36m888  888\e[0m \e[1;36m888 d88P\e[0m  |  DOI: 10.1016/j.softx.2024.101700
\e[1;36m888\e[0m        \e[1;36m"Y8888\e[0m  \e[1;36m888\e[0m     \e[1;36m888\e[0m \e[1;36m88888888\e[0m \e[1;36m"Y888888\e[0m \e[1;36m88888P"\e[0m   |  GitHub: https://github.com/PeriHub/PeriLab.jl\n""",
    )
end

"""
    parse_commandline()

Parse command-line arguments using the ArgParse package.

This function sets up argument parsing options for various command-line arguments and returns the parsed arguments as a dictionary.

## Arguments

None.

## Returns

- `Dict{String, Any}`: A dictionary containing the parsed command-line arguments.

## Command-Line Arguments

- `--dry_run`: If provided, it stores `true` in the dictionary.
- `--output_dir` or `-o`: If provided, it stores `true` in the dictionary.
- `--verbose` or `-v`: If provided, it stores `true` in the dictionary.
- `--debug` or `-d`: If provided, it stores `true` in the dictionary.
- `--silent` or `-s`: If provided, it stores `true` in the dictionary.
- `--reload` or `-r`: If provided, it stores `true` in the dictionary.
- `filename`: A required argument that stores the provided filename in the dictionary.

## Usage

```julia
parsed_args = parse_commandline()
"""

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--dry_run"
        help = "dry_run"
        action = :store_true
        "--output_dir", "-o"
        help = "output_dir"
        arg_type = String
        default = ""
        "--verbose", "-v"
        help = "verbose"
        action = :store_true
        "--debug", "-d"
        help = "debug"
        action = :store_true
        "--silent", "-s"
        help = "silent"
        action = :store_true
        "--reload", "-r"
        help = "reload"
        action = :store_true
        "filename"
        help = "filename"
        required = true
    end

    return parse_args(s)
end

"""
    main()

Entry point for the PeriLab application.

This function serves as the entry point for the PeriLab application. It calls the core `main` function with the provided arguments.
"""
function main()::Cint
    parsed_args = parse_commandline()
    @debug "Parsed args:"
    for (arg, val) in parsed_args
        @debug "  $arg  =>  $val"
    end
    MPI.Init()
    main(
        parsed_args["filename"];
        output_dir = parsed_args["output_dir"],
        dry_run = parsed_args["dry_run"],
        verbose = parsed_args["verbose"],
        debug = parsed_args["debug"],
        silent = parsed_args["silent"],
        reload = parsed_args["reload"],
    )
    return 0
end

"""
    get_examples()

Copy the examples folder to the current directory.
"""
function get_examples()
    package_dir = dirname(dirname(pathof(PeriLab)))
    examples_dir = joinpath(package_dir, "examples")
    dest_folder = joinpath(pwd(), "examples")
    if isdir(examples_dir)
        cp(examples_dir, dest_folder)
    end
end

"""
    main(filename::String, output_dir::String="", dry_run::Bool=false, verbose::Bool=false, debug::Bool=false, silent::Bool=false, reload::Bool=false)

Entry point for the PeriLab application.

This function serves as the entry point for the PeriLab application. It calls the core `main` function with the provided arguments.

# Arguments
- `filename::String`: The filename of the input file.
- `output_dir::String`: The output directory.
- `dry_run::Bool=false`: Whether to run in dry-run mode.
- `verbose::Bool=false`: Whether to run in verbose mode.
- `debug::Bool=false`: Whether to run in debug mode.
- `silent::Bool=false`: Whether to run in silent mode.
- `reload::Bool=false`: Whether to reload the input file.
"""
function main(
    filename::String;
    output_dir::String = "",
    dry_run::Bool = false,
    verbose::Bool = false,
    debug::Bool = false,
    silent::Bool = false,
    reload::Bool = false,
)

    @timeit to "PeriLab" begin
        if !MPI.Initialized()
            MPI.Init()
        end
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)

        result_files = nothing
        outputs = nothing

        try
            # atexit(() -> cleanup(comm))
            if silent && debug
                @warn "Silent and debug mode currently cannot be used at the same time."
                silent = false
            end
            Logging_module.init_logging(filename, debug, silent, rank, size)
            if rank == 0
                if !silent
                    print_banner()
                end
                @info "\n PeriLab version: $PERILAB_VERSION\n Copyright: Dr.-Ing. Christian Willberg, M.Sc. Jan-Timo Hesse\n Contact: christian.willberg@dlr.de, jan-timo.hesse@dlr.de\n GitHub: https://github.com/PeriHub/PeriLab.jl\n DOI: 10.1016/j.softx.2024.101700\n License: BSD-3-Clause\n ---------------------------------------------------------------\n"
                @info Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
                try
                    dirty, git_info =
                        Logging_module.get_current_git_info(joinpath(@__DIR__, ".."))
                    if dirty
                        @warn git_info
                    else
                        @info git_info
                    end
                catch e
                    if !isa(e, LibGit2.GitError)
                        @warn "No current git info."
                    end
                end
                if size > 1
                    @info "MPI: Running on " * string(size) * " processes"
                end
            end

            ################################
            filedirectory = dirname(filename)
            if output_dir != ""
                if !isdir(output_dir)
                    try
                        mkpath(output_dir)
                    catch
                        @error "Could not create output directory"
                    end
                end
            else
                output_dir = filedirectory
            end

            if !reload
                Data_manager.initialize_data()
            else
                @info "PeriLab started in the reload mode"
            end
            Data_manager.set_silent(silent)
            @timeit to "IO.initialize_data" datamanager, params =
                IO.initialize_data(filename, filedirectory, Data_manager, comm, to)
            @info "Init solver"
            @timeit to "Solver.init" block_nodes, bcs, datamanager, solver_options =
                Solver.init(params, datamanager, to)
            IO.show_block_summary(
                solver_options,
                params,
                Logging_module.get_log_file(),
                silent,
                comm,
                datamanager,
            )
            IO.show_mpi_summary(Logging_module.get_log_file(), silent, comm, datamanager)
            @debug "Init write results"
            @timeit to "IO.init_write_results" result_files, outputs =
                IO.init_write_results(
                    params,
                    output_dir,
                    filedirectory,
                    datamanager,
                    solver_options["nsteps"],
                    PERILAB_VERSION,
                )
            Logging_module.set_result_files(result_files)
            if dry_run
                nsteps = solver_options["nsteps"]
                solver_options["nsteps"] = 10
                elapsed_time = @elapsed begin
                    @timeit to "Solver" result_files = Solver.solver(
                        solver_options,
                        block_nodes,
                        bcs,
                        datamanager,
                        outputs,
                        result_files,
                        IO.write_results,
                        to,
                        silent,
                    )
                end

                @info "Estimated runtime: " * string((elapsed_time / 10) * nsteps) * " [s]"
                file_size = IO.get_file_size(result_files)
                @info "Estimated filesize: " * string((file_size / 10) * nsteps) * " [b]"

            else
                @timeit to "Solver.solver" result_files = Solver.solver(
                    solver_options,
                    block_nodes,
                    bcs,
                    datamanager,
                    outputs,
                    result_files,
                    IO.write_results,
                    to,
                    silent,
                )
            end

        catch e
            if e isa InterruptException
                @info "PeriLab was interrupted"
            elseif !isa(e, DomainError)
                rethrow(e)
            end
        end
        if !isnothing(result_files)
            @debug "Close result files"
            IO.close_result_files(result_files, outputs)

            if size > 1 && rank == 0
                IO.merge_exodus_files(result_files, output_dir)
            end
            MPI.Barrier(comm)
            if size > 1 || dry_run
                IO.delete_files(result_files, output_dir)
            end
        end
    end
    if verbose
        TimerOutputs.complement!(to)
        @info to
    end
    @info "PeriLab finished"
end

# PrecompileTools.@compile_workload begin
#     main("examples/Small/Input.yaml"; silent=true)
# end

end # module
