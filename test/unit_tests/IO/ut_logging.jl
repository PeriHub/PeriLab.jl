# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
# include("../../../src/IO/logging.jl")
#using Test
using LoggingExtras
# import .Logging_Module

@testset "ut_init_logging" begin
    PeriLab.Logging_Module.init_logging("test", false, false, 0, 1)
    # @test typeof(current_logger()) ==
    #       LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
    #                                                                         typeof(PeriLab.Logging_Module.progress_filter)},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#11#17"},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#12#18"},
    #                                                    Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    PeriLab.Logging_Module.init_logging("test", false, false, 0, 2)
    # @test typeof(current_logger()) ==
    #       LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
    #                                                                         typeof(PeriLab.Logging_Module.progress_filter)},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#11#17"},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#12#18"},
    #                                                    Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    PeriLab.Logging_Module.init_logging("test", false, false, 1, 2)
    # @test typeof(current_logger()) ==
    #       LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
    #                                                                         typeof(PeriLab.Logging_Module.progress_filter)},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#11#17"},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#12#18"},
    #                                                    Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    PeriLab.Logging_Module.init_logging("test", true, false, 0, 1)
    # @test typeof(current_logger()) ==
    #       LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
    #                                                                         typeof(PeriLab.Logging_Module.progress_filter)},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#7#13"},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#8#14"},
    #                                                    Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    PeriLab.Logging_Module.init_logging("test", true, false, 0, 2)
    # @test typeof(current_logger()) ==
    #       LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
    #                                                                         typeof(PeriLab.Logging_Module.progress_filter)},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#7#13"},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#8#14"},
    #                                                    Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test_2.0")
    PeriLab.Logging_Module.init_logging("test", true, false, 1, 2)
    # @test typeof(current_logger()) ==
    #       LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
    #                                                                         typeof(PeriLab.Logging_Module.progress_filter)},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#7#13"},
    #                                                    Base.CoreLogging.LogLevel},
    #                                     MinLevelLogger{FormatLogger{PeriLab.Logging_Module.var"#8#14"},
    #                                                    Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test_2.1")
    PeriLab.Logging_Module.init_logging("test", false, false, 0, 1)
end
