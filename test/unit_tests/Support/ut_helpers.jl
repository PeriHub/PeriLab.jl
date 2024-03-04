# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/helpers.jl")
include("../../../src/Support/data_manager.jl")
using Reexport
@reexport using .Helpers
using ProgressBars
@testset "ut_interpolation" begin
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = [-1.0, 0.0, 7.0, 26.0, 63.0]  # x.^3 - 1.
    values_dict = Dict()
    values_dict["value"] = Helpers.interpolation(x, y)
    @test isapprox(Helpers.interpol_data([1.5, 2.5], values_dict["value"])[1], 2.375)
    @test isapprox(Helpers.interpol_data([1.5, 2.5], values_dict["value"])[2], 14.625)
    @test isapprox(Helpers.interpol_data(1.5, values_dict["value"]), 2.375)
    @test Helpers.interpol_data(-1, values_dict["value"]) == minimum(y)
    @test Helpers.interpol_data([-1, -8], values_dict["value"]) == [minimum(y), minimum(y)]
    @test Helpers.interpol_data(5, values_dict["value"]) == maximum(y)
end
@testset "ut_find_indices" begin
    @test Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 1) == [1, 2, 9]
    @test Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 2) == [3]
    @test Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 3) == [4, 5]
    @test Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 4) == [6, 7, 8]
    @test Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 5) == []
end

# Define a test case for the find_active function
@testset "ut_find_active" begin
    # Test case 1: Empty input
    @test isempty(Helpers.find_active(Bool[]))
    # Test case 2: All elements are active
    @test Helpers.find_active([true, true, true]) == [1, 2, 3]
    # Test case 3: No elements are active
    @test isempty(Helpers.find_active([false, false, false]))
    # Test case 4: Mix of active and inactive elements
    @test Helpers.find_active([false, true, false, true, true]) == [2, 4, 5]
    # Test case 5: SubList of active and inactive elements
    list = [false, true, false, true, true]
    @test Helpers.find_active(list[[2, 3, 5]]) == [1, 3]
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
    inverse_nlist = Helpers.find_inverse_bond_id(nlist)

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
    @test Helpers.check_inf_or_nan(a, "a") == false
    a[1, 1] = 1 / 0
    @test Helpers.check_inf_or_nan(a, "Testing infinite test vector")
    a = 0
    @test Helpers.check_inf_or_nan(a, "a") == false
end
@testset "get_matrix_style" begin
    A = 1
    @test length(size(A)) == 0
    Atest = Helpers.matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1]
    @test length(size(A)) == 1
    @test sum(size(A)) == 1
    Atest = Helpers.matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1 1; 1 1]
    @test length(size(A)) == 2
    @test sum(size(A)) == 4
    Atest = Helpers.matrix_style(A)
    @test length(size(A)) == 2
    @test sum(size(A)) == 4
    A = [1 1 1; 1 1 1; 1 1 1]
    @test sum(size(A)) == 6
    Atest = Helpers.matrix_style(A)
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
    @test Helpers.find_files_with_ending(tmpdir, ".txt") == ["file1.txt", "file2.txt"]

    # Test case 2: Find .csv files
    @test Helpers.find_files_with_ending(tmpdir, ".csv") == ["file3.csv", "file4.csv"]

    # Clean up: Remove the temporary test directory and files
    rm(tmpdir; recursive=true)
end

# only interface test, because the called fromVoigt function is tested in "Tensors"
@testset "ut_get_fourth_order" begin
    @test size(Helpers.get_fourth_order(zeros(Float64, 6, 6), 3)) == (3, 3, 3, 3)
    @test size(Helpers.get_fourth_order(zeros(Float64, 3, 3), 2)) == (2, 2, 2, 2)
end

@testset "ut_progress_bar" begin
    nsteps::Int64 = rand(1:100)
    @test Helpers.progress_bar(rand(1:100), nsteps, true) == 1:nsteps+1
    @test Helpers.progress_bar(rand(1:100), nsteps, false) == 1:nsteps+1
    @test Helpers.progress_bar(0, nsteps, true) == 1:nsteps+1
    @test typeof(Helpers.progress_bar(0, nsteps, false)) == ProgressBar
    @test length(Helpers.progress_bar(0, nsteps, false)) == nsteps + 1
end

@testset "get_active_update_nodes" begin
    nnodes = 4
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(nnodes)
    update_list = test_Data_manager.create_constant_node_field("Update List", Bool, 1, true)
    active = test_Data_manager.create_constant_node_field("Active List", Bool, 1, true)
    block_nodes = Dict(1 => [1, 2], 2 => [3, 4])
    block = 1
    @test Helpers.get_active_update_nodes(active, update_list, block_nodes, block) == ([1, 2], [1, 2],)
    block = 2
    @test Helpers.get_active_update_nodes(active, update_list, block_nodes, block) == ([3, 4], [3, 4])
    update_list[3] = false
    @test Helpers.get_active_update_nodes(active, update_list, block_nodes, block) == ([3, 4], [4])
    active[3] = false
    @test Helpers.get_active_update_nodes(active, update_list, block_nodes, block) == ([4], [4])
    update_list[3] = true
    @test Helpers.get_active_update_nodes(active, update_list, block_nodes, block) == ([4], [4])
    update_list[3] = false
    update_list[4] = false
    @test Helpers.get_active_update_nodes(active, update_list, block_nodes, block) == ([4], [])
    active[3] = false
    active[4] = false
    @test Helpers.get_active_update_nodes(active, update_list, block_nodes, block) == ([], [])
end