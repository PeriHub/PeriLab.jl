# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/helpers.jl")
include("../../../src/Support/data_manager.jl")
using .Helpers
using ProgressBars
@testset "ut_interpolation" begin
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = [-1.0, 0.0, 7.0, 26.0, 63.0]  # x.^3 - 1.
    values_dict = Dict()
    values_dict["value"] = interpolation(x, y)
    @test interpol_data([1.5, 2.5], values_dict["value"]) == [2.375, 14.625]
    @test interpol_data(1.5, values_dict["value"]) == 2.375
    @test interpol_data(-1, values_dict["value"]) == minimum(y)
    @test interpol_data([-1, -8], values_dict["value"]) == [minimum(y), minimum(y)]
    @test interpol_data(5, values_dict["value"]) == maximum(y)
end
@testset "ut_find_indices" begin
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 1) == [1, 2, 9]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 2) == [3]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 3) == [4, 5]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 4) == [6, 7, 8]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 5) == []
end

# Define a test case for the find_active function
@testset "ut_find_active" begin
    # Test case 1: Empty input
    @test isempty(find_active(Bool[]))
    # Test case 2: All elements are active
    @test find_active([true, true, true]) == [1, 2, 3]
    # Test case 3: No elements are active
    @test isempty(find_active([false, false, false]))
    # Test case 4: Mix of active and inactive elements
    @test find_active([false, true, false, true, true]) == [2, 4, 5]
end

@testset "ut_find_inverse_bond_id" begin
    test_Data_manager = Data_manager
    nnodes = 2
    num_responder = 1
    test_Data_manager.set_num_controller(nnodes)
    test_Data_manager.set_num_responder(num_responder)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 3
    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1, 2]
    inverse_nlist = find_inverse_bond_id(nlist)

    @test length(inverse_nlist[1]) == 1
    @test length(inverse_nlist[2]) == 2
    @test length(inverse_nlist[3]) == 1
    @test inverse_nlist[1][2] == 1
    @test inverse_nlist[2][1] == 1
    @test inverse_nlist[2][3] == 2
    @test inverse_nlist[3][2] == 2

end

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

# only interface test, because the called fromVoigt function is tested in "Tensors"
@testset "ut_get_fourth_order" begin
    @test size(get_fourth_order(zeros(Float64, 6, 6), 3)) == (3, 3, 3, 3)
    @test size(get_fourth_order(zeros(Float64, 3, 3), 2)) == (2, 2, 2, 2)
end

@testset "ut_progress_bar" begin
    nsteps::Int64 = rand(1:100)
    @test progress_bar(rand(1:100), nsteps, true) == 1:nsteps+1
    @test progress_bar(rand(1:100), nsteps, false) == 1:nsteps+1
    @test progress_bar(0, nsteps, true) == 1:nsteps+1
    @test typeof(progress_bar(0, nsteps, false)) == ProgressBar
    @test length(progress_bar(0, nsteps, false)) == nsteps + 1
end