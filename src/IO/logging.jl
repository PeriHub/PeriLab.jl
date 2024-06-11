# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Logging_module
using Logging
using LoggingExtras
using TimerOutputs
using DataFrames
using LibGit2
using PrettyTables
using Dates
export init_logging
export set_result_files
export get_current_git_info

result_files::Vector{Dict} = []
log_file::String = ""

"""
    set_log_file(filename::String)

Set the log file.

# Arguments
- `filename::String`: The filename.
# Returns
- `log_file::String`: The log file.
"""
function set_log_file(filename::String, debug::Bool, rank::Int64, size::Int64)

    if size > 1
        if debug
            return split(filename, ".")[1] * "_" * Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS") * "_$size.$rank.log"
        elseif rank == 0
            return split(filename, ".")[1] * "_" * Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS") * ".log"
        else
            return ""
        end
    end
    return split(filename, ".")[1] * "_" * Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS") * ".log"

end

"""
    get_log_file()

Get the log file.

# Returns
- `log_file::String`: The log file.
"""
function get_log_file()
    global log_file
    return log_file
end

"""
    set_result_files(result_files_temp::Vector{Dict})

Set the result files.

# Arguments
- `result_files_temp::Vector{Dict}`: The result files.
"""
function set_result_files(result_files_temp::Vector{Dict})
    global result_files = result_files_temp
end

# function progress_wrap(f::Function, desc::String)
#     p = ProgressUnknown(desc=desc, spinner=true, color=:blue)
#     channel = Channel() do ch
#         while !isready(ch)
#             next!(p)
#         end
#         finish!(p)
#         take!(ch)
#     end

#     f()

#     put!(channel, 1)
# end

function close_result_file(result_file::Dict)
    if !isnothing(result_file["file"])
        close(result_file["file"])
    end
end

function close_result_files(result_files::Vector{Dict})
    for result_file in result_files
        try
            close_result_file(result_file)
            return true
        catch
            @warn "File already closed"
            return false
        end
    end
end

"""
    print_table(data::Matrix, datamanager::Module)

Print the table.

# Arguments
- `data::Matrix`: The data.
- `datamanager::Module`: The data manager.
"""
function print_table(data::Matrix, datamanager::Module)
    if !datamanager.get_silent()
        pretty_table(
            data;
            body_hlines=[1],
            body_hlines_format=Tuple('─' for _ = 1:4),
            cell_alignment=Dict((1, 1) => :l),
            formatters=ft_printf("%10.1f", 2),
            highlighters=(
                hl_cell([(1, 1)], crayon"bold"),
                hl_col(2, crayon"dark_gray")
            ),
            show_header=false,
            tf=tf_borderless
        )
        pretty_table(
            current_logger().loggers[2].logger.stream,
            data;
            body_hlines=[1],
            body_hlines_format=Tuple('─' for _ = 1:4),
            cell_alignment=Dict((1, 1) => :l),
            formatters=ft_printf("%10.1f", 2),
            highlighters=(
                hl_cell([(1, 1)], crayon"bold"),
                hl_col(2, crayon"dark_gray")
            ),
            show_header=false,
            tf=tf_borderless
        )
    else
        pretty_table(
            current_logger().loggers[1].logger.stream,
            data;
            body_hlines=[1],
            body_hlines_format=Tuple('─' for _ = 1:4),
            cell_alignment=Dict((1, 1) => :l),
            formatters=ft_printf("%10.1f", 2),
            highlighters=(
                hl_cell([(1, 1)], crayon"bold"),
                hl_col(2, crayon"dark_gray")
            ),
            show_header=false,
            tf=tf_borderless
        )
    end
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
    init_logging(filename::String, debug::Bool, silent::Bool, rank::Int64, size::Int64)

Initialize the logging.

# Arguments
- `filename::String`: The filename.
- `debug::Bool`: If debug is true.
- `silent::Bool`: If silent is true.
- `rank::Int64`: The rank.
- `size::Int64`: The size.
"""
function init_logging(filename::String, debug::Bool, silent::Bool, rank::Int64, size::Int64)

    global log_file

    log_file = set_log_file(filename, debug, rank, size)

    if log_file == ""
        Logging.disable_logging(Logging.Error)
        return
    end
    demux_logger = nothing
    try
        if debug
            file_logger = FormatLogger(log_file; append=false) do io, args
                if args.level in [Logging.Info, Logging.Warn, Logging.Error, Logging.Debug]
                    println(io, "[", args.level, "] ", args._module, ", ", args.line, " | ", args.message)
                end
            end
            filtered_logger = ActiveFilteredLogger(progress_filter, ConsoleLogger(stderr))
            demux_logger = TeeLogger(
                MinLevelLogger(filtered_logger, Logging.Debug),
                MinLevelLogger(file_logger, Logging.Debug),
            )
        elseif silent
            io = open(log_file, "a")
            redirect_stderr(io)
            file_logger = FormatLogger(log_file; append=false) do io, args
                if args.level in [Logging.Info, Logging.Warn, Logging.Error, Logging.Debug]
                    println(io, "[", args.level, "] ", args.message)
                end
            end
            error_logger = FormatLogger(log_file; append=false) do io, args
                if args.level == Logging.Error
                    close_result_files(result_files)
                    exit()
                end
            end
            demux_logger = TeeLogger(
                MinLevelLogger(file_logger, Logging.Debug),
                MinLevelLogger(error_logger, Logging.Info),
            )
        else
            file_logger = FormatLogger(log_file; append=false) do io, args
                if args.level in [Logging.Info, Logging.Warn, Logging.Error, Logging.Debug]
                    println(io, "[", args.level, "] ", args.message)
                end
            end
            error_logger = FormatLogger(log_file; append=false) do io, args
                if args.level == Logging.Error
                    close_result_files(result_files)
                    throw(DomainError(args, args.message))
                end
            end
            filtered_logger = ActiveFilteredLogger(progress_filter, ConsoleLogger(stderr))
            demux_logger = TeeLogger(
                MinLevelLogger(filtered_logger, Logging.Info),
                MinLevelLogger(file_logger, Logging.Debug),
                MinLevelLogger(error_logger, Logging.Info),
            )
        end
    catch e
        if e isa SystemError
            @error "Could not open log file: $log_file, make sure the directory exists."
            throw(DomainError(e))
        else
            rethrow(e)
        end
    end
    global_logger(demux_logger)
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