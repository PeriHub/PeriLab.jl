# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using LoggingExtras
import .Logging_module

#@testset "ut_init_logging" begin
Logging_module.init_logging("test", false, false, 0, 1)
@test typeof(current_logger()) == LoggingExtras.TeeLogger{Tuple{
    MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,typeof(Main.Logging_module.progress_filter)},Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel}}}
# @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
Logging_module.init_logging("test", false, false, 0, 2)
@test typeof(current_logger()) == LoggingExtras.TeeLogger{Tuple{
    MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,typeof(Main.Logging_module.progress_filter)},Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel}}}
# @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
Logging_module.init_logging("test", false, false, 1, 2)
@test typeof(current_logger()) == LoggingExtras.TeeLogger{Tuple{
    MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,typeof(Main.Logging_module.progress_filter)},Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel}}}
# @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
Logging_module.init_logging("test", true, false, 0, 1)
@test typeof(current_logger()) == LoggingExtras.TeeLogger{Tuple{
    MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,typeof(Main.Logging_module.progress_filter)},Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel}}}
# @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
Logging_module.init_logging("test", true, false, 0, 2)
@test typeof(current_logger()) == LoggingExtras.TeeLogger{Tuple{
    MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,typeof(Main.Logging_module.progress_filter)},Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel}}}
# @test contains(current_logger().loggers[2].logger.stream.name, "<file test_2.0")
Logging_module.init_logging("test", true, false, 1, 2)
@test typeof(current_logger()) == LoggingExtras.TeeLogger{Tuple{
    MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,typeof(Main.Logging_module.progress_filter)},Base.CoreLogging.LogLevel},
    MinLevelLogger{FormatLogger,Base.CoreLogging.LogLevel}}}
# @test contains(current_logger().loggers[2].logger.stream.name, "<file test_2.1")
#end