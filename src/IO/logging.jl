# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Logging_module
using Logging
using LoggingExtras
using TimerOutputs
using DataFrames
include("./IO.jl")
export init_logging
export set_result_files

result_files = []

function set_result_files(result_files_temp)
    global result_files = result_files_temp
end

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

function init_logging(filename, debug, rank, size)

    logfilename = split(filename, ".")[1] * ".log"

    if debug
        if size > 1
            logfilename = split(filename, ".")[1] * "$size.$rank.log"
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
end