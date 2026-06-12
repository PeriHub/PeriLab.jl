# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module PeriLabExceptions
using Logging
"""
    PeriLabError <: Exception

A custom exception used to stop execution silently when caught
by the main handler, while ensuring the error message is still logged.
"""
struct PeriLabError <: Exception
    val::Any
    msg::AbstractString
    PeriLabError(@nospecialize(val)) = (@noinline; new(val, ""))
    PeriLabError(@nospecialize(val), @nospecialize(msg)) = (@noinline; new(val, msg))
end

"""
    @silent_error(msg)

Logs the error message to your configured loggers and then stops
execution by throwing a `PeriLabError`.
"""
macro abort(msg)
    ex = esc(msg) # Evaluate in caller's scope, safely inject literal value
    quote
        Logging.@error($ex)
        throw(PeriLabError("Fatal", $ex))
    end
end
end
