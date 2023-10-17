# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/IO/mesh_data.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
using .Read_Mesh
using .Data_manager
using DataFrames

@testset "ut_bondIntersectInfinitePlane_2d" begin
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
        intersect_inf_plane, x = Read_Mesh.bondIntersectInfinitePlane(data[:, 1], data[:, i], lower_left_corner, normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end

    intersect_inf_plane, x = Read_Mesh.bondIntersectInfinitePlane(data[:, 6], data[:, 5], lower_left_corner, normal)
    @test intersect_inf_plane == false
    @test x == undef

    lower_left_corner = [0.0, 0.0]
    normal = [0.0, -1.0]
    #first value not important


    for i in 2:5
        intersect_inf_plane, x = Read_Mesh.bondIntersectInfinitePlane(data[:, 1], data[:, i], lower_left_corner, normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end
    lower_left_corner = [10.0, 0.0]
    normal = [0.0, 1.0]
    #first value not important
    for i in 2:5
        intersect_inf_plane, x = Read_Mesh.bondIntersectInfinitePlane(data[:, 1], data[:, i], lower_left_corner, normal)
        @test intersect_inf_plane == test_vals[i]
        @test x == test_coor[i]
    end
end


#bondIntersect