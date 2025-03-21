# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/IO/logging.jl")
using Test
using LoggingExtras
import .Logging_module

@testset "ut_init_logging" begin
    Logging_module.init_logging("test", false, false, 0, 1)
    @test typeof(current_logger()) ==
          LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
                                                                            typeof(Main.Logging_module.progress_filter)},
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    Logging_module.init_logging("test", false, false, 0, 2)
    @test typeof(current_logger()) ==
          LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
                                                                            typeof(Main.Logging_module.progress_filter)},
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    Logging_module.init_logging("test", false, false, 1, 2)
    @test typeof(current_logger()) ==
          LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
                                                                            typeof(Main.Logging_module.progress_filter)},
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    Logging_module.init_logging("test", true, false, 0, 1)
    @test typeof(current_logger()) ==
          LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
                                                                            typeof(Main.Logging_module.progress_filter)},
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test")
    Logging_module.init_logging("test", true, false, 0, 2)
    @test typeof(current_logger()) ==
          LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
                                                                            typeof(Main.Logging_module.progress_filter)},
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test_2.0")
    Logging_module.init_logging("test", true, false, 1, 2)
    @test typeof(current_logger()) ==
          LoggingExtras.TeeLogger{Tuple{MinLevelLogger{ActiveFilteredLogger{ConsoleLogger,
                                                                            typeof(Main.Logging_module.progress_filter)},
                                                       Base.CoreLogging.LogLevel},
                                        MinLevelLogger{FormatLogger,
                                                       Base.CoreLogging.LogLevel}}}
    # @test contains(current_logger().loggers[2].logger.stream.name, "<file test_2.1")
end

@testset "close_result_file tests" begin
    # Create a sample dictionary with an open file
    file_path = "test_output.txt"
    file = open(file_path, "w")
    result_file = Dict("file" => file)

    # Write something to the file
    write(result_file["file"], "Test data")

    # Call the function to close the file
    check = Logging_module.close_result_file(result_file)
    @test check
    # Check if the file is closed
    is_closed = !isopen(result_file["file"])

    # Test if the file is closed
    @test is_closed

    # Clean up
    rm(file_path, force = true)
    # Create a sample dictionary with file set to nothing
    result_file = Dict("file" => nothing)

    # Call the function to close the file
    check = Logging_module.close_result_file(result_file)
    @test !check
    # Test if the file is still nothing
    @test isnothing(result_file["file"])
end

@testset "close_result_files tests" begin
    # Create sample dictionaries with open files
    file_path1 = "test_output1.txt"
    file1 = open(file_path1, "w")
    result_file1 = Dict("file" => file1)

    file_path2 = "test_output2.txt"
    file2 = open(file_path2, "w")
    result_file2 = Dict("file" => file2)

    result_files = [result_file1, result_file2]

    # Write something to the files
    write(result_file1["file"], "Test data 1")
    write(result_file2["file"], "Test data 2")

    # Call the function to close the files
    result = Logging_module.close_result_files(result_files)

    # Check if the files are closed
    is_closed1 = !isopen(result_file1["file"])
    is_closed2 = !isopen(result_file2["file"])

    # Test if the files are closed and the function returns true
    @test is_closed1 == true
    @test is_closed2 == true
    @test result == true

    # Clean up
    rm(file_path1, force = true)
    rm(file_path2, force = true)
end
