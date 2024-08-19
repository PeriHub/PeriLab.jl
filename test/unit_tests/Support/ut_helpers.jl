# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using ProgressBars

#include("../../../src/PeriLab.jl")
#using .PeriLab

@testset "ut_invert" begin
    @test isnothing(PeriLab.Solver.Helpers.invert(zeros(Float64, 3, 3)))
    @test isnothing(PeriLab.Solver.Helpers.invert(zeros(Int64, 3, 3)))
    @test isnothing(
        PeriLab.Solver.Helpers.invert(zeros(Float64, 3, 3), "test Float is singular."),
    )
    @test isnothing(
        PeriLab.Solver.Helpers.invert(zeros(Int64, 3, 3), "test Int is singular."),
    )
    @test PeriLab.Solver.Helpers.invert([1 0; 0 1]) == inv([1 0; 0 1])
    @test PeriLab.Solver.Helpers.invert(
        [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
        "test Int is singular.",
    ) == inv([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
end


@testset "ut_find_local_neighbors" begin
    coordinates = zeros(Float64, 5, 2)
    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[3, 1] = 0
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1
    coordinates[4, 2] = 1
    coordinates[5, 1] = 2
    coordinates[5, 2] = 2

    nlist = [2, 3, 4, 5]

    bond_horizon::Float64 = 1

    @test PeriLab.Solver.Helpers.find_local_neighbors(
        5,
        coordinates,
        nlist,
        bond_horizon,
    ) == []
    bond_horizon = 2.6
    @test PeriLab.Solver.Helpers.find_local_neighbors(
        5,
        coordinates,
        nlist,
        bond_horizon,
    ) == [2, 3, 4]
    @test PeriLab.Solver.Helpers.find_local_neighbors(
        2,
        coordinates,
        nlist,
        bond_horizon,
    ) == [3, 4, 5]
end
@testset "ut_qdim" begin
    @test isnothing(PeriLab.Solver.Helpers.qdim(0, 2))
    @test isnothing(PeriLab.Solver.Helpers.qdim(0, 3))
    @test PeriLab.Solver.Helpers.qdim(1, 2) == 2
    @test PeriLab.Solver.Helpers.qdim(2, 2) == 5
    @test PeriLab.Solver.Helpers.qdim(3, 2) == 9
    @test PeriLab.Solver.Helpers.qdim(4, 2) == 14
    @test PeriLab.Solver.Helpers.qdim(5, 2) == 20
    @test PeriLab.Solver.Helpers.qdim(6, 2) == 27
    @test PeriLab.Solver.Helpers.qdim(7, 2) == 35
    @test PeriLab.Solver.Helpers.qdim(8, 2) == 44
    @test PeriLab.Solver.Helpers.qdim(9, 2) == 54
    @test PeriLab.Solver.Helpers.qdim(10, 2) == 65

    @test PeriLab.Solver.Helpers.qdim(1, 3) == 3
    @test PeriLab.Solver.Helpers.qdim(2, 3) == 9
    @test PeriLab.Solver.Helpers.qdim(3, 3) == 19
    @test PeriLab.Solver.Helpers.qdim(4, 3) == 34
    @test PeriLab.Solver.Helpers.qdim(5, 3) == 55
    @test PeriLab.Solver.Helpers.qdim(6, 3) == 83
    @test PeriLab.Solver.Helpers.qdim(7, 3) == 119
    @test PeriLab.Solver.Helpers.qdim(8, 3) == 164
    @test PeriLab.Solver.Helpers.qdim(9, 3) == 219
    @test PeriLab.Solver.Helpers.qdim(10, 3) == 285
end


@testset "rotate_second_order_tensor" begin

    angles = [0]
    rot = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensor = zeros(2, 2)
    tensor[1, 1] = 1
    dof = 2
    back = true
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        back,
    )
    @test tensorTest == tensor
    angles = [90.0]
    rot = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        back,
    )
    @test isapprox(tensorTest[1, 1] + 1, 1) # plus one, because of how approx works
    @test isapprox(tensorTest[1, 2] + 1, 1)
    @test isapprox(tensorTest[2, 1] + 1, 1)
    # @test isapprox(tensorTest[2, 2], 1)
    back = false
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        back,
    )
    @test tensorTest == tensor

    angles = [0, 0, 0]
    rot = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensor = zeros(3, 3)
    tensor[1, 1] = 1
    dof = 3

    back = true
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        back,
    )
    @test tensorTest == tensor
    angles = [0, 0, 90.0]
    rot = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        back,
    )
    @test isapprox(tensorTest[1, 1] + 1, 1) # plus one, because of how approx works
    @test isapprox(tensorTest[1, 2] + 1, 1)
    @test isapprox(tensorTest[1, 3] + 1, 1)
    @test isapprox(tensorTest[2, 1] + 1, 1)
    # @test isapprox(tensorTest[2, 2], 1)
    @test isapprox(tensorTest[2, 3] + 1, 1)
    @test isapprox(tensorTest[3, 1] + 1, 1)
    @test isapprox(tensorTest[3, 2] + 1, 1)
    @test isapprox(tensorTest[3, 3] + 1, 1)

    back = false
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        back,
    )
    @test tensorTest == tensor

    angles = [10, 20, 90.0]
    rot = PeriLab.IO.Geometry.rotation_tensor(angles)
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        true,
    )
    tensorTest = PeriLab.Solver.Helpers.rotate_second_order_tensor(
        Matrix{Float64}(rot),
        tensor,
        false,
    )
    @test tensorTest == tensor

end

@testset "ut_interpolation" begin
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = [-1.0, 0.0, 7.0, 26.0, 63.0]  # x.^3 - 1.
    values_dict = Dict()
    values_dict["value"] = PeriLab.Solver.Helpers.interpolation(x, y)
    @test isapprox(
        PeriLab.Solver.Helpers.interpol_data([1.5, 2.5], values_dict["value"])[1],
        2.375,
    )
    @test isapprox(
        PeriLab.Solver.Helpers.interpol_data([1.5, 2.5], values_dict["value"])[2],
        14.625,
    )
    @test isapprox(PeriLab.Solver.Helpers.interpol_data(1.5, values_dict["value"]), 2.375)
    @test PeriLab.Solver.Helpers.interpol_data(-1, values_dict["value"]) == minimum(y)
    @test PeriLab.Solver.Helpers.interpol_data([-1, -8], values_dict["value"]) ==
          [minimum(y), minimum(y)]
    @test PeriLab.Solver.Helpers.interpol_data(5, values_dict["value"]) == maximum(y)
    x = [0.0, 1.0]
    y = [-1.0, 0.0]
    values_dict = Dict()
    values_dict["value"] = PeriLab.Solver.Helpers.interpolation(x, y)
end
@testset "ut_find_indices" begin
    @test PeriLab.Solver.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 1) == [1, 2, 9]
    @test PeriLab.Solver.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 2) == [3]
    @test PeriLab.Solver.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 3) == [4, 5]
    @test PeriLab.Solver.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 4) == [6, 7, 8]
    @test PeriLab.Solver.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 5) == []
end

# Define a test case for the find_active function
@testset "ut_find_active" begin
    # Test case 1: Empty input
    @test isempty(PeriLab.Solver.Helpers.find_active(Bool[]))
    # Test case 2: All elements are active
    @test PeriLab.Solver.Helpers.find_active([true, true, true]) == [1, 2, 3]
    # Test case 3: No elements are active
    @test isempty(PeriLab.Solver.Helpers.find_active([false, false, false]))
    # Test case 4: Mix of active and inactive elements
    @test PeriLab.Solver.Helpers.find_active([false, true, false, true, true]) == [2, 4, 5]
    # Test case 5: SubList of active and inactive elements
    list = [false, true, false, true, true]
    @test PeriLab.Solver.Helpers.find_active(list[[2, 3, 5]]) == [1, 3]
end

@testset "ut_find_inverse_bond_id" begin
    test_data_manager = PeriLab.Data_manager
    nnodes = 2
    num_responder = 1
    test_data_manager.set_num_controller(nnodes)
    test_data_manager.set_num_responder(num_responder)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 3
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1, 2]
    inverse_nlist = PeriLab.Solver.Helpers.find_inverse_bond_id(nlist)

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
    @test PeriLab.Solver.Helpers.check_inf_or_nan(a, "a") == false
    a[1, 1] = 1 / 0
    @test PeriLab.Solver.Helpers.check_inf_or_nan(a, "Testing infinite test vector")
    a = 0
    @test PeriLab.Solver.Helpers.check_inf_or_nan(a, "a") == false
    a = NaN
    @test PeriLab.Solver.Helpers.check_inf_or_nan(a, "a") == true
end
@testset "get_matrix_style" begin
    A = 1
    @test length(size(A)) == 0
    Atest = PeriLab.Solver.Helpers.matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1]
    @test length(size(A)) == 1
    @test sum(size(A)) == 1
    Atest = PeriLab.Solver.Helpers.matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1 1; 1 1]
    @test length(size(A)) == 2
    @test sum(size(A)) == 4
    Atest = PeriLab.Solver.Helpers.matrix_style(A)
    @test length(size(A)) == 2
    @test sum(size(A)) == 4
    A = [1 1 1; 1 1 1; 1 1 1]
    @test sum(size(A)) == 6
    Atest = PeriLab.Solver.Helpers.matrix_style(A)
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
    @test PeriLab.Solver.Helpers.find_files_with_ending(tmpdir, ".txt") ==
          ["file1.txt", "file2.txt"]

    # Test case 2: Find .csv files
    @test PeriLab.Solver.Helpers.find_files_with_ending(tmpdir, ".csv") ==
          ["file3.csv", "file4.csv"]

    # Clean up: Remove the temporary test directory and files
    rm(tmpdir; recursive = true)
end

# only interface test, because the called fromVoigt function is tested in "Tensors"
@testset "ut_get_fourth_order" begin
    @test size(PeriLab.Solver.Helpers.get_fourth_order(zeros(Float64, 6, 6), 3)) ==
          (3, 3, 3, 3)
    @test size(PeriLab.Solver.Helpers.get_fourth_order(zeros(Float64, 3, 3), 2)) ==
          (2, 2, 2, 2)
end

@testset "ut_progress_bar" begin
    nsteps::Int64 = rand(1:100)
    @test PeriLab.Solver.Helpers.progress_bar(rand(1:100), nsteps, true) == 1:nsteps+1
    @test PeriLab.Solver.Helpers.progress_bar(rand(1:100), nsteps, false) == 1:nsteps+1
    @test PeriLab.Solver.Helpers.progress_bar(0, nsteps, true) == 1:nsteps+1
    @test typeof(PeriLab.Solver.Helpers.progress_bar(0, nsteps, false)) == ProgressBar
    @test length(PeriLab.Solver.Helpers.progress_bar(0, nsteps, false)) == nsteps + 1
end

@testset "get_active_update_nodes" begin
    nnodes = 4
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_num_controller(nnodes)
    update_list = test_data_manager.create_constant_node_field("Update List", Bool, 1, true)
    active = test_data_manager.create_constant_node_field("Active List", Bool, 1, true)
    block_nodes = Dict(1 => [1, 2], 2 => [3, 4])
    block = 1
    @test PeriLab.Solver.Helpers.get_active_update_nodes(
        active,
        update_list,
        block_nodes,
        block,
    ) == ([1, 2], [1, 2])
    block = 2
    @test PeriLab.Solver.Helpers.get_active_update_nodes(
        active,
        update_list,
        block_nodes,
        block,
    ) == ([3, 4], [3, 4])
    update_list[3] = false
    @test PeriLab.Solver.Helpers.get_active_update_nodes(
        active,
        update_list,
        block_nodes,
        block,
    ) == ([3, 4], [4])
    active[3] = false
    @test PeriLab.Solver.Helpers.get_active_update_nodes(
        active,
        update_list,
        block_nodes,
        block,
    ) == ([4], [4])
    update_list[3] = true
    @test PeriLab.Solver.Helpers.get_active_update_nodes(
        active,
        update_list,
        block_nodes,
        block,
    ) == ([4], [4])
    update_list[3] = false
    update_list[4] = false
    @test PeriLab.Solver.Helpers.get_active_update_nodes(
        active,
        update_list,
        block_nodes,
        block,
    ) == ([4], [])
    active[3] = false
    active[4] = false
    @test PeriLab.Solver.Helpers.get_active_update_nodes(
        active,
        update_list,
        block_nodes,
        block,
    ) == ([], [])
end
