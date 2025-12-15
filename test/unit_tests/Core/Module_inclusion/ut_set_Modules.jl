# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module MockModule
struct TestType end
get_name() = "existing_function"
call_function(args...) = sum(args)
end

using Test
using Random
include("../../../../src/Core/Module_inclusion/set_Modules.jl")

@testset "ut_find_jl_files" begin
    Random.seed!(rand(1:100000))
    base = "test_tmp_Set_modules"

    @test isnothing(find_jl_files(base))

    if isdir(base)
        rm(base, recursive = true)
    end

    mkdir(base)
    folder = base * "/" * randstring(12)
    mkpath(folder)
    subfolder1 = folder * "/" * randstring(12)
    mkpath(subfolder1)
    subfolder2 = folder * "/" * randstring(2)
    mkpath(subfolder2)
    subsubfolder = subfolder2 * "/" * randstring(3)
    mkpath(subsubfolder)
    filename1 = randstring(3)
    filename2 = randstring(8)
    filename2 = randstring(6)
    filename3 = randstring(2)
    filename4 = randstring(1)
    filename5 = randstring(8)
    io = open(folder * "/" * filename1 * ".jl", "w")
    close(io)
    io = open(subfolder1 * "/" * filename2 * ".jl", "w")
    close(io)
    io = open(subfolder1 * "/" * filename3 * ".jl", "w")
    close(io)
    io = open(subfolder1 * "/" * filename4 * ".jl", "w")
    close(io)
    io = open(subsubfolder * "/" * filename5 * ".jl", "w")
    close(io)
    io = open(folder * "/" * filename4 * ".dat", "w")
    close(io)

    list = find_jl_files(base)
    folder * "/" * filename1 * ".jl" in list
    @test subfolder1 * "/" * filename2 * ".jl" in list
    @test subfolder1 * "/" * filename3 * ".jl" in list
    @test subfolder1 * "/" * filename4 * ".jl" in list
    @test subsubfolder * "/" * filename5 * ".jl" in list
    @test !(folder * "/" * filename4 * ".dat" in list)
    rm(base, recursive = true)
end

@testset "ut_create_module_specifics" begin

    # Mock module setup

    filename::AbstractString = "MockModule.jl"

    module_list = [
        Dict{String,AbstractString}("File" => filename,
                                    "Module Name" => "MockModule"),
    ]

    specifics = Dict("Name" => "get_name",
                     "Call Function" => "call_function")

    # This should trigger the error because "nonexistent_function"
    # doesn't match "existing_function"
    result = create_module_specifics("nonexistent_function",
                                     module_list,
                                     @__MODULE__,
                                     specifics,
                                     (1, 2, 3))

    @test result === nothing

    result = create_module_specifics("nonexistent_function",
                                     module_list,
                                     @__MODULE__,
                                     "get_name")

    @test result === nothing
end
function create_test_file(filepath::String, content::String)
    open(filepath, "w") do f
        write(f, content)
    end
end
@testset "ut_find_module_files" begin
    mktempdir() do dir
        # Create a test file with module and function
        test_file = joinpath(dir, "TestModule.jl")
        content = """
        module TestModule

        export my_function

        function my_function()
        	return "test"
        end

        end
        """
        create_test_file(test_file, content)

        # Test finding the function
        result = find_module_files(dir, "my_function")

        @test length(result) == 1
        @test result[1]["File"] == test_file
        @test result[1]["Module Name"] == "TestModule"
    end
end
@testset "ut_find_module_files_not_found" begin
    mktempdir() do dir
        # Create a test file without the target function
        test_file = joinpath(dir, "TestModule.jl")
        content = """
        module TestModule

        function other_function()
        	return "test"
        end

        end
        """
        create_test_file(test_file, content)

        # Test - should return empty list
        result = find_module_files(dir, "my_function")

        @test length(result) == 0
    end
end
