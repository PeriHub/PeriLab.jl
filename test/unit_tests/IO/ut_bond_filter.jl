# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#using DataFrames

#using PeriLab

@testset "ut_bond_intersect_infinite_plane_2d" begin
    data = zeros(Float64, 2, 6)
    data[1, 1] = 0
    data[2, 1] = -0.5

    data[1, 2] = 0
    data[2, 2] = 0.5

    data[1, 3] = 0
    data[2, 3] = -1.5

    data[1, 4] = 1.0
    data[2, 4] = -0.5

    data[1, 5] = 0.5
    data[2, 5] = 0.0

    data[1, 6] = 1.0
    data[2, 6] = 0.0
    lower_left_corner = [0.0, 0.0]
    normal = [0.0, 1.0]

    test_vals = [true, true, false, false, true]
    test_coor = [undef, [0.0, 0.0], undef, undef, [0.5, 0.0]]
    #first value not important
    for i in 2:5
        intersect_inf_plane,
        x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 1],
                                                     data[:, i],
                                                     lower_left_corner,
                                                     normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end

    intersect_inf_plane,
    x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 6],
                                                 data[:, 5],
                                                 lower_left_corner,
                                                 normal)
    @test intersect_inf_plane == false
    @test x == undef

    lower_left_corner = [0.0, 0.0]
    normal = [0.0, -1.0]
    #first value not important

    for i in 2:5
        intersect_inf_plane,
        x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 1],
                                                     data[:, i],
                                                     lower_left_corner,
                                                     normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end
    lower_left_corner = [10.0, 0.0]
    normal = [0.0, 1.0]
    #first value not important
    for i in 2:5
        intersect_inf_plane,
        x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 1],
                                                     data[:, i],
                                                     lower_left_corner,
                                                     normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end
end

@testset "ut_bond_intersect_infinite_plane_3d" begin
    data = zeros(Float64, 3, 6)
    data[1, 1] = 0
    data[2, 1] = -0.5

    data[1, 2] = 0
    data[2, 2] = 0.5

    data[1, 3] = 0
    data[2, 3] = -1.5

    data[1, 4] = 1.0
    data[2, 4] = -0.5

    data[1, 5] = 0.5
    data[2, 5] = 0.0

    data[1, 6] = 1.0
    data[2, 6] = 0.0
    lower_left_corner = [0.0, 0.0, 0.0]
    normal = [0.0, 1.0, 0.0]

    test_vals = [true, true, false, false, true]
    test_coor = [undef, [0.0, 0.0, 0.0], undef, undef, [0.5, 0.0, 0.0]]
    #first value not important
    for i in 2:5
        intersect_inf_plane,
        x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 1],
                                                     data[:, i],
                                                     lower_left_corner,
                                                     normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end

    intersect_inf_plane,
    x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 6],
                                                 data[:, 5],
                                                 lower_left_corner,
                                                 normal)
    @test intersect_inf_plane == false
    @test x == undef

    lower_left_corner = [0.0, 0.0, 0.0]
    normal = [0.0, -1.0, 0.0]
    #first value not important
    for i in 2:5
        intersect_inf_plane,
        x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 1],
                                                     data[:, i],
                                                     lower_left_corner,
                                                     normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end
    lower_left_corner = [10.0, 0.0, 5.0]
    normal = [0.0, 1.0, 0.0]
    #first value not important
    for i in 2:5
        intersect_inf_plane,
        x = PeriLab.IO.bond_intersect_infinite_plane(data[:, 1],
                                                     data[:, i],
                                                     lower_left_corner,
                                                     normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end
end
# y is sorted out, because all points which are not in the infinite plane are not included
@testset "ut_bondIntersect" begin
    lower_left_corner = [0.0, 0.0, 0.0]
    bottom_unit_vector = [1.0, 0.0, 0.0]
    normal = [0.0, 1.0, 0.0]
    side_length = 1.0
    bottom_length = 1.0
    x = [0.0, 1.0, 0.0]
    bond_intersect = PeriLab.IO.bond_intersect_rectangle_plane(x,
                                                               lower_left_corner,
                                                               bottom_unit_vector,
                                                               normal,
                                                               side_length,
                                                               bottom_length)
    @test bond_intersect == true
    x = [0.0, 0.0, 0.0]
    bond_intersect = PeriLab.IO.bond_intersect_rectangle_plane(x,
                                                               lower_left_corner,
                                                               bottom_unit_vector,
                                                               normal,
                                                               side_length,
                                                               bottom_length)
    @test bond_intersect == true
    x = [10.0, 0.0, 0.0]
    bond_intersect = PeriLab.IO.bond_intersect_rectangle_plane(x,
                                                               lower_left_corner,
                                                               bottom_unit_vector,
                                                               normal,
                                                               side_length,
                                                               bottom_length)
    @test bond_intersect == false
    x = [0.0, 0.0, 5.0]
    bond_intersect = PeriLab.IO.bond_intersect_rectangle_plane(x,
                                                               lower_left_corner,
                                                               bottom_unit_vector,
                                                               normal,
                                                               side_length,
                                                               bottom_length)
    @test bond_intersect == false
    x = [-0.2, 0.0, 0.0]
    bond_intersect = PeriLab.IO.bond_intersect_rectangle_plane(x,
                                                               lower_left_corner,
                                                               bottom_unit_vector,
                                                               normal,
                                                               side_length,
                                                               bottom_length)
    @test bond_intersect == false
    normal = [0.0, -1.0, 0.0]
    bond_intersect = PeriLab.IO.bond_intersect_rectangle_plane(x,
                                                               lower_left_corner,
                                                               bottom_unit_vector,
                                                               normal,
                                                               side_length,
                                                               bottom_length)
    @test bond_intersect == false
    normal = [0.0, -1.0, 0.0]
    bottom_unit_vector = [-1.0, 0.0, 0.0]
    bond_intersect = PeriLab.IO.bond_intersect_rectangle_plane(x,
                                                               lower_left_corner,
                                                               bottom_unit_vector,
                                                               normal,
                                                               side_length,
                                                               bottom_length)
    @test bond_intersect == true
end

@testset "ut_disk_filter" begin
    nnodes = 8

    data = zeros(Float64, 3, 8)
    data[1, 1] = 0
    data[2, 1] = 0
    data[3, 1] = -1.0

    data[1, 2] = 0
    data[2, 2] = 0.5
    data[3, 2] = -1.0

    data[1, 3] = 0
    data[2, 3] = -0.5
    data[3, 3] = -1.0

    data[1, 4] = 2.0
    data[2, 4] = -0.5
    data[3, 4] = -1.0

    data[1, 5] = 0
    data[2, 5] = 0
    data[3, 5] = 1.0

    data[1, 6] = 0
    data[2, 6] = 0.5
    data[3, 6] = 1.0

    data[1, 7] = 0
    data[2, 7] = -0.5
    data[3, 7] = 1.0

    data[1, 8] = 2.0
    data[2, 8] = -0.5
    data[3, 8] = 1.0

    filter = Dict("Center X" => 0.0,
                  "Center Y" => 0.0,
                  "Center Z" => 0.0,
                  "Normal X" => 0.0,
                  "Normal Y" => 0.0,
                  "Normal Z" => 1.0,
                  "Radius" => 1.0)

    nlist = [
        [2, 3, 4, 5, 6, 7, 8],
        [1, 3, 4, 5, 6, 7, 8],
        [1, 2, 4, 5, 6, 7, 8],
        [1, 2, 3, 5, 6, 7, 8],
        [1, 2, 3, 4, 6, 7, 8],
        [1, 2, 3, 4, 5, 7, 8],
        [1, 2, 3, 4, 5, 6, 8],
        [1, 2, 3, 4, 5, 6, 7]
    ]
    dof = 3

    # Define the expected output values
    expected_filter_flag = [
        [true, true, true, false, false, false, true],   # Node 1
        [true, true, true, false, false, false, true],   # Node 2
        [true, true, true, false, false, false, true],   # Node 3
        [true, true, true, true, true, true, true],   # Node 4
        [false, false, false, true, true, true, true],   # Node 5
        [false, false, false, true, true, true, true],    # Node 6
        [false, false, false, true, true, true, true],   # Node 7
        [true, true, true, true, true, true, true]    # Node 8
    ]

    expected_normal = [0, 0, 1]
    (filter_flag, normal) = PeriLab.IO.disk_filter(nnodes, data, filter, nlist, dof)
    @test filter_flag == expected_filter_flag
    @test normal == expected_normal

    @test isnothing(PeriLab.IO.disk_filter(nnodes, data, filter, nlist, 2))
end

@testset "ut_rectangular_plane_filter" begin
    nnodes = 6

    data = zeros(Float64, 3, 6)
    data[1, 1] = 0
    data[2, 1] = 0
    data[3, 1] = -1.0

    data[1, 2] = 0
    data[2, 2] = -0.5
    data[3, 2] = -1.0

    data[1, 3] = 0
    data[2, 3] = 1.5
    data[3, 3] = -1.0

    data[1, 4] = 0
    data[2, 4] = 0
    data[3, 4] = 1.0

    data[1, 5] = 0
    data[2, 5] = -0.5
    data[3, 5] = 1.0

    data[1, 6] = 0
    data[2, 6] = 1.5
    data[3, 6] = 1.0

    filter = Dict("Lower Left Corner X" => -0.5,
                  "Lower Left Corner Y" => -0.5,
                  "Lower Left Corner Z" => 0.0,
                  "Bottom Unit Vector X" => 1.0,
                  "Bottom Unit Vector Y" => 0.0,
                  "Bottom Unit Vector Z" => 0.0,
                  "Normal X" => 0.0,
                  "Normal Y" => 0.0,
                  "Normal Z" => 1.0,
                  "Bottom Length" => 1.0,
                  "Side Length" => 1.0)

    nlist = [
        [2, 3, 4, 5, 6],
        [1, 3, 4, 5, 6],
        [1, 2, 4, 5, 6],
        [1, 2, 3, 5, 6],
        [1, 2, 3, 4, 6],
        [1, 2, 3, 4, 5]
    ]
    dof = 3

    # Define the expected output values
    expected_filter_flag = [
        [true, true, true, true, true],   # Node 1
        [true, true, true, false, true],   # Node 2
        [true, true, true, true, true],   # Node 3
        [true, true, true, true, true],   # Node 4
        [true, false, true, true, true],   # Node 5
        [true, true, true, true, true]    # Node 6
    ]

    expected_normal = [0, 0, 1]
    (filter_flag,
     normal) = PeriLab.IO.rectangular_plane_filter(nnodes, data, filter, nlist,
                                                   dof)
    @test filter_flag == expected_filter_flag
    @test normal == expected_normal
end
