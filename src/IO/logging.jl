# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Logging_module
using Logging
using LoggingExtras
using TimerOutputs
using DataFrames
using LibGit2
include("./IO.jl")
export init_logging
export set_result_files
export get_current_git_info

result_files::Vector{Dict} = []

"""
    set_result_files(result_files_temp::Vector{Dict})

Set the result files.

# Arguments
- `result_files_temp::Vector{Dict}`: The result files.
"""
function set_result_files(result_files_temp::Vector{Dict})
    global result_files = result_files_temp
end

"""
    progress_filter(log_args)

Filter progress messages.

# Arguments
- `log_args`: The log arguments.
# Returns
- `true`: If the message is not a progress message.
- `false`: If the message is a progress message.
"""
function progress_filter(log_args)
    if log_args.message isa TimerOutputs.TimerOutput || log_args.message isa DataFrames.DataFrame
        return true
    end
    if log_args.message isa String
        if startswith(log_args.message, "Step:")
            return false
        end
        if startswith(log_args.message, "\n PeriLab version:")
            return false
        end
    end
    return true
end

"""
    init_logging(filename::String, debug::Bool, rank::Int64, size::Int64)

Initialize the logging.

# Arguments
- `filename::String`: The filename.
- `debug::Bool`: If debug is true.
- `rank::Int64`: The rank.
- `size::Int64`: The size.
"""
function init_logging(filename::String, debug::Bool, rank::Int64, size::Int64)

    logfilename = split(filename, ".")[1] * ".log"

    if debug
        if size > 1
            logfilename = split(filename, ".")[1] * "_$size.$rank.log"
        end
        file_logger = FormatLogger(logfilename; append=false) do io, args
            if args.level in [Logging.Info, Logging.Warn, Logging.Error, Logging.Debug]
                println(io, "[", args.level, "] ", args._module, ", ", args.line, " | ", args.message)
            end
        end
        filtered_logger = ActiveFilteredLogger(progress_filter, ConsoleLogger(stderr))
        demux_logger = TeeLogger(
            MinLevelLogger(filtered_logger, Logging.Debug),
            MinLevelLogger(file_logger, Logging.Debug),
        )
        global_logger(demux_logger)
    else
        file_logger = FormatLogger(logfilename; append=false) do io, args
            if args.level in [Logging.Info, Logging.Warn, Logging.Error, Logging.Debug]
                println(io, "[", args.level, "] ", args.message)
            end
        end
        error_logger = FormatLogger(logfilename; append=false) do io, args
            if args.level == Logging.Error
                IO.close_result_files(result_files)
                exit()
            end
        end
        filtered_logger = ActiveFilteredLogger(progress_filter, ConsoleLogger(stderr))
        demux_logger = TeeLogger(
            MinLevelLogger(filtered_logger, Logging.Info),
            MinLevelLogger(file_logger, Logging.Info),
            MinLevelLogger(error_logger, Logging.Info),
        )

        if rank == 0
            global_logger(demux_logger)
        else
            Logging.disable_logging(Logging.Error)
        end
    end
end

function get_current_git_info(repo_path::AbstractString=".")
    repo = LibGit2.GitRepo(repo_path)
    head_name = LibGit2.headname(repo)
    head_oid = LibGit2.head_oid(repo)
    dirty = LibGit2.isdirty(repo)
    tag_list = LibGit2.tag_list(repo)

    info = ""

    if head_name == "main"
        info *= "Commit: $head_oid"
    else
        info *= "Branch: $head_name, Commit: $head_oid"
    end

    for tag in tag_list
        tag_hash = LibGit2.target(LibGit2.GitObject(repo, tag))
        if tag_hash == head_oid
            info *= ", Tag: $tag"
            break
        end
    end

    if dirty
        info *= ", Local changes detected"
    end
    return dirty, info
end
end