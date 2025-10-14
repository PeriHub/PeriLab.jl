# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
	PeriLab

A module for managing and executing peridynamic simulations in PeriLab.

This module provides functionality for running simulations in the PeriLab environment. It includes functions for initializing simulations, processing data, running solvers, and managing results.

## Modules

- `Core/Data_manager.jl`: Data manager module for data management and access.
- `IO/logging.jl`: Module for setting up and managing logging.
- `IO/IO.jl`: Input/output functions for handling data files.
- `Core/Solver/Solver_manager.jl`: Solver control module for managing simulation solvers.

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
include("./Support/Helpers.jl")
include("./Support/Geometry.jl")
include("./Core/Data_manager.jl")
include("./IO/logging.jl")
include("./MPI_communication/MPI_communication.jl")
include("./Support/Parameters/parameter_handling.jl")
include("./IO/IO.jl")
include("./Core/Solver/Solver_manager.jl")

using MPI
using TimerOutputs
using Logging
using ArgParse
using Dates
using LibGit2
using StyledStrings

const to = TimerOutput()

using .Data_Manager

import .Logging_Module
import .IO
using .Solver_Manager

PERILAB_VERSION = "1.5.0"

export main

"""
	print_banner(mpi)

Prints a banner displaying information about the PeriLab application.

This function prints a banner containing details about the PeriLab application, including its name, version, copyright, contact information, and license. It provides a visual introduction to the application.
"""
function print_banner(mpi)
    line = []
    push!(line,
          styled"{bright_blue:PeriLab.                   }{bright_green:d8b} {bright_blue:888               888       }| {bold:Version}: {underline:$PERILAB_VERSION}")
    push!(line,
          styled"{bright_blue:888   Y88b                 }{bright_green:Y8P} {bright_blue:888               888       }| {bold:Copyright}:")
    push!(line,
          styled"{bright_blue:888    888                     888               888       }|  Dr.-Ing. Christian Willberg (https://orcid.org/0000-0003-2433-9183)")
    push!(line,
          styled"{bright_blue:888   d88P .d88b.  888d888 888 888       8888b.  88888b.   }|  M.Sc. Jan-Timo Hesse (https://orcid.org/0000-0002-3006-1520)")
    push!(line,
          styled"{bright_blue:8888888P\" d8P  Y8b 888P\"   888 888          \"88b 888 \"88b  }| {bold:Contact}: christian.willberg@h2.de, jan-timo.hesse@dlr.de")
    push!(line,
          styled"{bright_blue:888       88888888 888     888 888      .d888888 888  888  }| {bold:GitHub}: https://github.com/PeriHub/PeriLab.jl")
    push!(line,
          styled"{bright_blue:888       Y8b.     888     888 888      888  888 888 d88P  }| {bold:DOI}: 10.1016/j.softx.2024.101700")
    push!(line,
          styled"{bright_blue:888        \"Y8888  888     888 88888888 \"Y888888 88888P\"   }| {bold:License}: BSD-3-Clause")

    num_of_chars = 0
    for l in line
        num_of_chars = max(num_of_chars, length(l) - displaysize(stdout)[2])
    end
    for l in line
        if num_of_chars == 0 || mpi
            println(l)
        else
            println(split(l, "|")[1][1:(end - num_of_chars)] * "|" * split(l, "|")[2])
        end
    end
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
        "filenames"
        nargs = '*'
        help = "filenames"
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
    for filename in parsed_args["filenames"]
        main(filename;
             output_dir = parsed_args["output_dir"],
             dry_run = parsed_args["dry_run"],
             verbose = parsed_args["verbose"],
             debug = parsed_args["debug"],
             silent = parsed_args["silent"],
             reload = parsed_args["reload"])
    end
    MPI.Finalize()
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
function main(filename::String;
              output_dir::String = "",
              dry_run::Bool = false,
              verbose::Bool = false,
              debug::Bool = false,
              silent::Bool = false,
              reload::Bool = false)
    reset_timer!(to)
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
            Logging_Module.init_logging(filename, debug, silent, rank, size)
            if rank == 0
                if !silent
                    print_banner(size > 1)
                end
                @info "\n PeriLab version: $PERILAB_VERSION\n Copyright: Dr.-Ing. Christian Willberg, M.Sc. Jan-Timo Hesse\n Contact: christian.willberg@dlr.de, jan-timo.hesse@dlr.de\n GitHub: https://github.com/PeriHub/PeriLab.jl\n DOI: 10.1016/j.softx.2024.101700\n License: BSD-3-Clause\n ---------------------------------------------------------------\n"
                @info Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
                try
                    dirty,
                    git_info = Logging_Module.get_current_git_info(joinpath(@__DIR__,
                                                                            ".."))
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
                Data_Manager.initialize_data()
            else
                @info "PeriLab started in the reload mode"
            end
            Data_Manager.set_silent(silent)
            Data_Manager.set_verbose(verbose)
            @timeit to "IO.initialize_data" datamanager,
                                            params,
                                            steps=IO.initialize_data(filename,
                                                                     filedirectory,
                                                                     Data_Manager,
                                                                     comm, to)
            datamanager.set_max_step(steps[end])
            for step_id in steps
                if !isnothing(step_id)
                    @info "Step: " * string(step_id) * " of " * string(length(steps))
                end
                datamanager.set_cancel(false)
                datamanager.set_step(step_id)
                @info "Init Solver"
                @timeit to "Solver_Manager.init" block_nodes,
                                                 bcs,
                                                 datamanager,
                                                 solver_options=Solver_Manager.init(params,
                                                                                    datamanager,
                                                                                    to,
                                                                                    step_id)
                if datamanager.get_current_time() >= solver_options["Final Time"]
                    @info "Step " * string(step_id) * " skipped."
                    continue
                end
                if datamanager.get_current_time() < solver_options["Initial Time"]
                    @info "Initial time not reached. Skipping step " * string(step_id)
                    continue
                end
                @timeit to "IO.init orientations" datamanager=IO.init_orientations(datamanager)
                IO.show_block_summary(solver_options,
                                      params,
                                      Logging_Module.get_log_file(),
                                      silent,
                                      comm,
                                      datamanager)
                IO.show_mpi_summary(Logging_Module.get_log_file(),
                                    silent,
                                    comm,
                                    datamanager)
                @debug "Init write results"
                if isnothing(step_id) || step_id == 1
                    @timeit to "IO.init_write_results" result_files,
                                                       outputs=IO.init_write_results(params,
                                                                                     output_dir,
                                                                                     filedirectory,
                                                                                     datamanager,
                                                                                     PERILAB_VERSION)
                end
                IO.set_output_frequency(params,
                                        datamanager,
                                        solver_options["Number of Steps"],
                                        step_id)
                if verbose
                    fields = datamanager.get_all_field_keys()
                    @info "Found " * string(length(fields)) * " Fields"
                    @info fields
                end
                if dry_run
                    nsteps = solver_options["Number of Steps"]
                    solver_options["Number of Steps"] = 10
                    elapsed_time = @elapsed begin
                        @timeit to "Solver" result_files=Solver_Manager.solver(solver_options,
                                                                               block_nodes,
                                                                               bcs,
                                                                               datamanager,
                                                                               outputs,
                                                                               result_files,
                                                                               IO.write_results,
                                                                               to,
                                                                               silent)
                    end

                    @info "Estimated runtime: " *
                          string((elapsed_time / 10) * nsteps) *
                          " [s]"
                    file_size = IO.get_file_size(result_files)
                    @info "Estimated filesize: " *
                          string((file_size / 10) * nsteps) *
                          " [b]"

                else
                    @timeit to "Solver_Manager.solver" result_files=Solver_Manager.solver(solver_options,
                                                                                          block_nodes,
                                                                                          bcs,
                                                                                          datamanager,
                                                                                          outputs,
                                                                                          result_files,
                                                                                          IO.write_results,
                                                                                          to,
                                                                                          silent)
                end
            end

        catch e
            if e isa InterruptException
                @info "PeriLab was interrupted"
            elseif !isa(e, Logging_Module.PeriLabError)
                open(Logging_Module.log_file, "a") do io
                    println(io, "[Error] ", e)
                end
                if !silent
                    rethrow(e)
                end
            end
        end
        if !isnothing(result_files)
            @debug "Close result files"
            IO.close_result_files(result_files, outputs)

            if size > 1 && rank == 0
                IO.merge_exodus_files(result_files, output_dir)
            end
            MPI.Barrier(comm)
            if (size > 1 && !debug) || dry_run
                IO.delete_files(result_files, output_dir)
            end
        end
    end
    if verbose
        TimerOutputs.complement!(to)
        @info to
        reset_timer!(to)
    end
    @info Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    @info "PeriLab finished"
end

end # module
