# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using ProgressBars
using LinearAlgebra
using StaticArrays
#using PeriLab

@testset "ut_create_permutation" begin
    nnodes = 3
    dof = 2
    perm = PeriLab.Solver_Manager.Helpers.create_permutation(nnodes, dof)
    @test perm == [1, 2, 3, 4, 5, 6]
    @test length(perm) == length(nodes) * dof
    nodes = [2, 5, 7]
    dof = 2
    perm = PeriLab.Solver_Manager.Helpers.create_permutation(nodes, dof, 7)
    @test perm == [2, 5, 7, 9, 12, 14]
    nodes = collect(1:100)
    dof = 3
    perm = PeriLab.Solver_Manager.Helpers.create_permutation(nodes, dof, 100)

    @test length(perm) == 300
    @test all(1 .<= perm .<= 300)

    perm = PeriLab.Solver_Manager.Helpers.create_permutation([1, 3], 2)
    @test length(perm) == 4
    @test perm == [1, 3, 4, 6]
end

@testset "ut_remove_ids" begin
    test_dict = Dict(1 => 1, 2 => 2, 3 => 4, 7 => 2)
    PeriLab.Solver_Manager.Helpers.remove_ids(test_dict, Vector{Int64}([1, 7, 8]))
    @test test_dict == Dict(2 => 2, 3 => 4)

    test_vec = [1, 2, 3, 7]
    PeriLab.Solver_Manager.Helpers.remove_ids(test_vec, Vector{Int64}([1, 7, 8]))
    @test test_vec == [2, 3]
end

@testset "ut_smat" begin
    @test PeriLab.Solver_Manager.Helpers.smat(zeros(2, 2)) == zeros(2, 2)
    @test PeriLab.Solver_Manager.Helpers.smat(zeros(3, 3)) == zeros(3, 3)
    @test PeriLab.Solver_Manager.Helpers.smat(zeros(4, 4)) == zeros(4, 4)
end

@testset "ut_determinant" begin
    @test PeriLab.Solver_Manager.Helpers.determinant(zeros(2, 2)) == 0
    @test PeriLab.Solver_Manager.Helpers.determinant(zeros(3, 3)) == 0
    @test PeriLab.Solver_Manager.Helpers.determinant(zeros(4, 4)) == 0
    @test PeriLab.Solver_Manager.Helpers.determinant(ones(3, 3)) == det(ones(3, 3))
end

@testset "ut_invert" begin
    test_inv = zeros(2, 2)
    PeriLab.Solver_Manager.Helpers.invert(zeros(Float64, 4, 4),
                                          "test Float is singular.")

    PeriLab.Solver_Manager.Helpers.invert(test_inv, [1.0 0; 0 1])
    @test test_inv == inv([1 0; 0 1])
    test_inv = zeros(3, 3)
    PeriLab.Solver_Manager.Helpers.invert(test_inv, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
                                          "test Int is singular.")
    @test test_inv == inv([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
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

    @test PeriLab.Solver_Manager.Helpers.find_local_neighbors(5,
                                                              coordinates,
                                                              nlist,
                                                              bond_horizon) == []
    bond_horizon = 2.6
    @test PeriLab.Solver_Manager.Helpers.find_local_neighbors(5,
                                                              coordinates,
                                                              nlist,
                                                              bond_horizon) == [2, 3, 4]
    @test PeriLab.Solver_Manager.Helpers.find_local_neighbors(2,
                                                              coordinates,
                                                              nlist,
                                                              bond_horizon) == [3, 4, 5]
end
@testset "ut_fastdot" begin
    a = rand(3)
    b = rand(3)
    PeriLab.Solver_Manager.Helpers.fastdot(a, b) == dot(a, b)
end
@testset "ut_get_mapping" begin
    @test length(PeriLab.Solver_Manager.Helpers.get_mapping(4)) == 0
    @test length(PeriLab.Solver_Manager.Helpers.get_mapping(1)) == 0
end

@testset "ut_qdim" begin
    @test isnothing(PeriLab.Solver_Manager.Helpers.qdim(0, 2))
    @test isnothing(PeriLab.Solver_Manager.Helpers.qdim(0, 3))
    @test PeriLab.Solver_Manager.Helpers.qdim(1, 2) == 2
    @test PeriLab.Solver_Manager.Helpers.qdim(2, 2) == 5
    @test PeriLab.Solver_Manager.Helpers.qdim(3, 2) == 9
    @test PeriLab.Solver_Manager.Helpers.qdim(4, 2) == 14
    @test PeriLab.Solver_Manager.Helpers.qdim(5, 2) == 20
    @test PeriLab.Solver_Manager.Helpers.qdim(6, 2) == 27
    @test PeriLab.Solver_Manager.Helpers.qdim(7, 2) == 35
    @test PeriLab.Solver_Manager.Helpers.qdim(8, 2) == 44
    @test PeriLab.Solver_Manager.Helpers.qdim(9, 2) == 54
    @test PeriLab.Solver_Manager.Helpers.qdim(10, 2) == 65

    @test PeriLab.Solver_Manager.Helpers.qdim(1, 3) == 3
    @test PeriLab.Solver_Manager.Helpers.qdim(2, 3) == 9
    @test PeriLab.Solver_Manager.Helpers.qdim(3, 3) == 19
    @test PeriLab.Solver_Manager.Helpers.qdim(4, 3) == 34
    @test PeriLab.Solver_Manager.Helpers.qdim(5, 3) == 55
    @test PeriLab.Solver_Manager.Helpers.qdim(6, 3) == 83
    @test PeriLab.Solver_Manager.Helpers.qdim(7, 3) == 119
    @test PeriLab.Solver_Manager.Helpers.qdim(8, 3) == 164
    @test PeriLab.Solver_Manager.Helpers.qdim(9, 3) == 219
    @test PeriLab.Solver_Manager.Helpers.qdim(10, 3) == 285
end

@testset "rotate_second_order_tensor" begin
    angles = [0]
    rot = PeriLab.Geometry.rotation_tensor(angles, 2)
    tensor = zeros(2, 2)
    tensor[1, 1] = 1
    dof = 2
    back = true
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           back)
    @test tensorTest == tensor
    angles = [90.0]
    rot = PeriLab.Geometry.rotation_tensor(angles, 2)
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           back)

    @test isapprox(tensorTest[1, 1] + 1, 1) # plus one, because of how approx works
    @test isapprox(tensorTest[1, 2] + 1, 1)
    @test isapprox(tensorTest[2, 1] + 1, 1)
    @test isapprox(tensorTest[2, 2] + 1, 1)
    back = false
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           back)
    @test tensorTest == tensor

    angles = [0, 0, 0]
    rot = PeriLab.Geometry.rotation_tensor(angles, 3)
    tensor = zeros(3, 3)
    tensor[1, 1] = 1
    dof = 3

    back = true
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           back)
    @test tensorTest == tensor
    angles = [0, 0, 90.0]
    rot = PeriLab.Geometry.rotation_tensor(angles, 3)
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           back)
    @test isapprox(tensorTest[1, 1] + 1, 1) # plus one, because of how approx works
    @test isapprox(tensorTest[1, 2] + 1, 1)
    @test isapprox(tensorTest[1, 3] + 1, 1)
    @test isapprox(tensorTest[2, 1] + 1, 1)
    @test isapprox(tensorTest[2, 2] + 1, 1)
    @test isapprox(tensorTest[2, 3] + 1, 1)
    @test isapprox(tensorTest[3, 1] + 1, 1)
    @test isapprox(tensorTest[3, 2] + 1, 1)
    @test isapprox(tensorTest[3, 3] + 1, 1)

    back = false
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           back)
    @test tensorTest == tensor

    angles = [10, 20, 90.0]
    rot = PeriLab.Geometry.rotation_tensor(angles, 3)
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           true)
    tensorTest = PeriLab.Solver_Manager.Helpers.rotate_second_order_tensor(Matrix{Float64}(rot),
                                                                           tensor,
                                                                           false)
    @test tensorTest == tensor
end

@testset "ut_interpolation" begin
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = [-1.0, 0.0, 7.0, 26.0, 63.0]  # x.^3 - 1.
    values_dict = Dict()
    values_dict["value"] = PeriLab.Solver_Manager.Helpers.interpolation(x, y)
    @test isapprox(PeriLab.Solver_Manager.Helpers.interpol_data([1.5, 2.5],
                                                                values_dict["value"])[1],
                   2.375)
    @test isapprox(PeriLab.Solver_Manager.Helpers.interpol_data([1.5, 2.5],
                                                                values_dict["value"])[2],
                   14.625)
    @test isapprox(PeriLab.Solver_Manager.Helpers.interpol_data(1.5, values_dict["value"]),
                   2.375)
    @test PeriLab.Solver_Manager.Helpers.interpol_data(-1, values_dict["value"]) ==
          minimum(y)
    @test PeriLab.Solver_Manager.Helpers.interpol_data([-1, -8], values_dict["value"]) ==
          [minimum(y), minimum(y)]
    @test PeriLab.Solver_Manager.Helpers.interpol_data(5, values_dict["value"]) ==
          maximum(y)
    x = [0.0, 1.0]
    y = [-1.0, 0.0]
    values_dict = Dict()
    values_dict["value"] = PeriLab.Solver_Manager.Helpers.interpolation(x, y)
end
@testset "ut_find_indices" begin
    @test PeriLab.Solver_Manager.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 1) ==
          [1, 2, 9]
    @test PeriLab.Solver_Manager.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 2) == [3]
    @test PeriLab.Solver_Manager.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 3) ==
          [4, 5]
    @test PeriLab.Solver_Manager.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 4) ==
          [6, 7, 8]
    @test PeriLab.Solver_Manager.Helpers.find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 5) == []
end

# Define a test case for the find_active_nodes function
@testset "ut_find_active_nodes" begin
    index = [0, 0, 0, 0, 0]
    # Test case 1: Empty input
    @test isempty(PeriLab.Solver_Manager.Helpers.find_active_nodes(Bool[], index, 1:0))
    # Test case 2: All elements are active
    @test PeriLab.Solver_Manager.Helpers.find_active_nodes([true, true, true],
                                                           index,
                                                           1:3) == [1, 2, 3]
    # Test case 3: No elements are active
    @test isempty(PeriLab.Solver_Manager.Helpers.find_active_nodes([false, false, false],
                                                                   index, 1:3))

    @test PeriLab.Solver_Manager.Helpers.find_active_nodes([false, false, false],
                                                           index,
                                                           1:3,
                                                           false) == [1, 2, 3]
    # Test case 4: Mix of active and inactive elements
    @test PeriLab.Solver_Manager.Helpers.find_active_nodes([false, true, false, true, true],
                                                           index,
                                                           1:5) == [2, 4, 5]
    # Test case 5: SubList of active and inactive elements
    list = [false, true, false, true, true]
    @test PeriLab.Solver_Manager.Helpers.find_active_nodes(list[[2, 3, 5]], index, 1:3) ==
          [1, 3]
end

@testset "ut_find_inverse_bond_id" begin
    test_data_manager = PeriLab.Data_Manager
    nnodes = 2
    num_responder = 1
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nnodes)
    test_data_manager.set_num_responder(num_responder)
    nn = test_data_manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 3
    nlist = test_data_manager.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1, 2]
    inverse_nlist = PeriLab.Solver_Manager.Helpers.find_inverse_bond_id(nlist)

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
    @test PeriLab.Solver_Manager.Helpers.check_inf_or_nan(a, "a") == false
    a[1, 1] = 1 / 0
    @test PeriLab.Solver_Manager.Helpers.check_inf_or_nan(a, "Testing infinite test vector")
    a[1, 1] = NaN
    @test PeriLab.Solver_Manager.Helpers.check_inf_or_nan(a, "Testing infinite test vector")
    a = 0
    @test PeriLab.Solver_Manager.Helpers.check_inf_or_nan(a, "a") == false
    a = NaN
    @test PeriLab.Solver_Manager.Helpers.check_inf_or_nan(a, "a")
    a = 1 / 0
    @test PeriLab.Solver_Manager.Helpers.check_inf_or_nan(a, "a") == true
end
@testset "get_matrix_style" begin
    A = 1
    Atest = PeriLab.Solver_Manager.Helpers.matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1;;]
    Atest = PeriLab.Solver_Manager.Helpers.matrix_style(A)
    @test sum(size(Atest)) == 2
    A = [1 1; 1 1]
    Atest = PeriLab.Solver_Manager.Helpers.matrix_style(A)
    @test length(size(Atest)) == 2
    @test sum(size(Atest)) == 4
    A = [1 1 1; 1 1 1; 1 1 1]
    Atest = PeriLab.Solver_Manager.Helpers.matrix_style(A)
    @test sum(size(Atest)) == 6
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
    @test PeriLab.Solver_Manager.Helpers.find_files_with_ending(tmpdir, ".txt") ==
          ["file1.txt", "file2.txt"]

    # Test case 2: Find .csv files
    @test PeriLab.Solver_Manager.Helpers.find_files_with_ending(tmpdir, ".csv") ==
          ["file3.csv", "file4.csv"]

    # Clean up: Remove the temporary test directory and files
    rm(tmpdir; recursive = true)
end

# only interface test, because the called fromVoigt function is tested in "Tensors"
@testset "ut_get_fourth_order" begin
    @test size(PeriLab.Solver_Manager.Helpers.get_fourth_order(zeros(Float64, 6, 6),
                                                               3)) == (3, 3, 3, 3)
    @test size(PeriLab.Solver_Manager.Helpers.get_fourth_order(zeros(Float64, 3, 3),
                                                               2)) == (2, 2, 2, 2)
    @test size(PeriLab.Solver_Manager.Helpers.get_fourth_order(zeros(Float64, 3, 3),
                                                               1)) == (0, 0, 0, 0)
end

@testset "ut_progress_bar" begin
    nsteps::Int64 = rand(1:100)
    @test PeriLab.Solver_Manager.Helpers.progress_bar(rand(1:100), nsteps, true) ==
          1:(nsteps + 1)
    @test PeriLab.Solver_Manager.Helpers.progress_bar(rand(1:100), nsteps, false) ==
          1:(nsteps + 1)
    @test PeriLab.Solver_Manager.Helpers.progress_bar(0, nsteps, true) == 1:(nsteps + 1)
    @test typeof(PeriLab.Solver_Manager.Helpers.progress_bar(0, nsteps, false)) ==
          ProgressBar
    @test length(PeriLab.Solver_Manager.Helpers.progress_bar(0, nsteps, false)) ==
          nsteps + 1
end

@testset "ut_get_active_update_nodes" begin
    nnodes = 4
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nnodes)
    update_list = test_data_manager.create_constant_node_scalar_field("Update", Bool;
                                                                      default_value = true)
    active = test_data_manager.create_constant_node_scalar_field("Active", Bool;
                                                                 default_value = true)
    update_nodes = test_data_manager.create_constant_node_scalar_field("Update Nodes",
                                                                       Int64)
    block_nodes = Dict(1 => [1, 2], 2 => [3, 4])
    block = 1
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [1, 2]
    block = 2
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [3, 4]
    update_list[3] = false
    block = 1
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [1, 2]
    block = 2
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [4]

    active[3] = false
    block = 1
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [1, 2]
    block = 2
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [4]
    update_list[3] = true
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [4]
    block = 1
    update_list[3] = false
    update_list[4] = false
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == [1, 2]
    block = 2
    active[3] = false
    active[4] = false
    @test PeriLab.Solver_Manager.Helpers.get_active_update_nodes(active,
                                                                 update_list,
                                                                 block_nodes[block],
                                                                 update_nodes) == []
end

@testset "ut_get_MMatrix" begin
    @test PeriLab.Solver_Manager.Helpers.get_MMatrix(4) ==
          MMatrix{2,2}(zeros(Float64, 2, 2))
    @test PeriLab.Solver_Manager.Helpers.get_MMatrix(9) ==
          MMatrix{3,3}(zeros(Float64, 3, 3))
    @test PeriLab.Solver_Manager.Helpers.get_MMatrix(36) ==
          MMatrix{6,6}(zeros(Float64, 6, 6))
    @test size(PeriLab.Solver_Manager.Helpers.get_MMatrix(8)) == (1, 1)
end

@testset "ut_sub_in_place" begin
    C = [[[0.0, 0], [0, 0]], [[0, 0], [0, 0]]]
    A = [[[1, 2], [3, 4]], [[5, 6.1], [7, 8]]]
    B = [[[1.0, 4], [5, 1]], [[2, 1], [1, 4]]]
    PeriLab.Solver_Manager.Helpers.sub_in_place!(C, A, B)
    @test C == [[[0, -2.0], [-2.0, 3.0]], [[3.0, 5.1], [6.0, 4.0]]]
end

@testset "ut_add_in_place" begin
    C = [[0.0, 0], [0, 0], [0, 0], [0, 0]]
    A = [[1, 2], [3, 4], [5, 6.1], [7, 8]]
    B = [[1.0, 4], [5, 1], [2, 1], [1, 4]]
    PeriLab.Solver_Manager.Helpers.add_in_place!(C, A, B)
    @test C == [[2.0, 6.0], [8.0, 5.0], [7.0, 7.1], [8.0, 12.0]]
end

@testset "ut_div_in_place" begin
    C = [[0.0, 0], [0, 0], [0, 0], [0, 0]]
    A = [[1, 2], [3, 4], [5, 6.1], [7, 8]]
    B = [[1.0, 4], [5, 1], [2, 1], [1, 4]]
    PeriLab.Solver_Manager.Helpers.div_in_place!(C, A, B)
    @test C == [[1.0, 0.5], [0.6, 4.0], [2.5, 6.1], [7.0, 2.0]]

    C = [0.0, 0.0]
    A = [-2.0, 4.0]
    B = 2.0
    PeriLab.Solver_Manager.Helpers.div_in_place!(C, A, B, false)
    @test C == [-1.0, 2.0]
    PeriLab.Solver_Manager.Helpers.div_in_place!(C, A, B, true)
    @test C == [1.0, 2.0]
end

@testset "ut_is_dependent" begin
    (fieldN, fieldNP1) = test_data_manager.create_node_scalar_field("Parameter", Float64)
    params = Dict("Value" => Dict("Field" => "Parameter"))
    @test PeriLab.Solver_Manager.Helpers.is_dependent("Value", params) ==
          (true, fieldNP1)
    params = Dict("Value" => 1.0)
    @test PeriLab.Solver_Manager.Helpers.is_dependent("Value", params) ==
          (false, nothing)
    params = Dict("Value" => Dict("Field" => "Non_Existent"))
    @test isnothing(PeriLab.Solver_Manager.Helpers.is_dependent("Value", params))
end

@testset "ut_matrix_to_voigt" begin
    matrix = Matrix{Float64}([1 2; 3 4])
    voigt = PeriLab.Solver_Manager.Helpers.matrix_to_voigt(matrix)
    @test voigt[1] == 1
    @test voigt[2] == 4
    @test voigt[3] == 2.5
    matrix = Matrix{Float64}([1 2 3; 4 5 6; 7 8 9])
    voigt = PeriLab.Solver_Manager.Helpers.matrix_to_voigt(matrix)
    @test voigt[1] == 1
    @test voigt[2] == 5
    @test voigt[3] == 9
    @test voigt[4] == 7
    @test voigt[5] == 5
    @test voigt[6] == 3
    matrix = Matrix{Float64}([1 2 3 3; 4 5 6 3; 7 8 9 3])
    @test isnothing(PeriLab.Solver_Manager.Helpers.matrix_to_voigt(matrix))
end
@testset "ut_voigt_to_matrix" begin
    @test isnothing(PeriLab.Solver_Manager.Helpers.voigt_to_matrix([1, 2.2]))
end
