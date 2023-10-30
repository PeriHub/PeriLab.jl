# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Support/tools.jl")

using Test

@testset "ut_check_inf_or_nan" begin
    a = ones(2, 2)
    @test check_inf_or_nan(a, "a") == false
    a[1, 1] = 1 / 0
    @test check_inf_or_nan(a, "Testing infinite test vector")
    a = 0
    @test check_inf_or_nan(a, "a") == false
end
@testset "get_matrix_style" begin
    A = 1
    @test length(size(A)) == 0
    Atest = matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1]
    @test length(size(A)) == 1
    @test sum(size(A)) == 1
    Atest = matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1 1; 1 1]
    @test length(size(A)) == 2
    @test sum(size(A)) == 4
    Atest = matrix_style(A)
    @test length(size(A)) == 2
    @test sum(size(A)) == 4
    A = [1 1 1; 1 1 1; 1 1 1]
    @test sum(size(A)) == 6
    Atest = matrix_style(A)
    @test sum(size(A)) == 6
end

@testset "ut_find_files_with_ending" begin
    # Create a temporary test directory with sample files
    tmpdir = mktempdir()
    touch(joinpath(tmpdir, "file1.txt"))
    touch(joinpath(tmpdir, "file2.txt"))
    touch(joinpath(tmpdir, "file3.csv"))
    touch(joinpath(tmpdir, "file4.csv"))
    mkdir(joinpath(tmpdir, "subdir"))


    # Test case 1: Find .txt files
    @test find_files_with_ending(tmpdir, ".txt") == ["file1.txt", "file2.txt"]

    # Test case 2: Find .csv files
    @test find_files_with_ending(tmpdir, ".csv") == ["file3.csv", "file4.csv"]

    # Clean up: Remove the temporary test directory and files
    rm(tmpdir; recursive=true)
end


@testset "ut_progress_bar" begin
    nsteps::Int64 = rand(1:100)
    @test progress_bar(rand(1:100), nsteps, true) == 1:nsteps+1
    @test progress_bar(rand(1:100), nsteps, false) == 1:nsteps+1
    @test progress_bar(0, nsteps, true) == 1:nsteps+1
    @test typeof(progress_bar(0, nsteps, false)) == ProgressBar
    @test length(progress_bar(0, nsteps, false)) == nsteps + 1
end

# only interface test, because the called fromVoigt function is tested in "Tensorial"
@testset "ut_get_fourth_order" begin
    @test size(get_fourth_order(zeros(Float64, 6, 6), 3)) == (3, 3, 3, 3)
    @test size(get_fourth_order(zeros(Float64, 3, 3), 2)) == (2, 2, 2, 2)
end

