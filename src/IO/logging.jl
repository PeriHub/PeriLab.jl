# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Logging_module
using Logging
using LoggingExtras
using TimerOutputs
export init_logging

function progress_filter(log_args)
    if typeof(log_args.message) == TimerOutputs.TimerOutput
        return true
    end
    if startswith(log_args.message, "Step:")
        return false
    end
    if startswith(log_args.message, "\n PeriLab version:")
        return false
    end
    return true
end

function init_logging(filename, debug)
    file_logger = FormatLogger(split(filename, ".")[1] * ".log"; append=false) do io, args
        if args.level in [Logging.Info, Logging.Warn, Logging.Error, Logging.Debug]
            if debug
                println(io, args._module, " | ", "[", args.level, "] ", args.message)
            else
                println(io, "[", args.level, "] ", args.message)
            end
        end
    end
    filtered_logger = ActiveFilteredLogger(progress_filter, ConsoleLogger(stderr))
    if debug
        demux_logger = TeeLogger(
            MinLevelLogger(file_logger, Logging.Debug),
            MinLevelLogger(filtered_logger, Logging.Debug),
        )
    else
        demux_logger = TeeLogger(
            MinLevelLogger(file_logger, Logging.Info),
            MinLevelLogger(filtered_logger, Logging.Info),
        )
    end
    global_logger(demux_logger)
end
end