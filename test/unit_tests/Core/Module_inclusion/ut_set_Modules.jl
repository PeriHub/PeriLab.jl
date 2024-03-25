# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


include("../../../../src/Core/Module_inclusion/set_Modules.jl")

using Test
using Random
# using .Set_modules

@testset "ut_find_jl_files" begin
    Random.seed!(rand(1:100000))
    base = "test_tmp_Set_modules"

    @test isnothing(Set_modules.find_jl_files(base))

    if isdir(base)
        rm(base, recursive=true)
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

    list = PeriLab.Solver.FEM.Set_modules.find_jl_files(base)
    folder * "/" * filename1 * ".jl" in list
    @test subfolder1 * "/" * filename2 * ".jl" in list
    @test subfolder1 * "/" * filename3 * ".jl" in list
    @test subfolder1 * "/" * filename4 * ".jl" in list
    @test subsubfolder * "/" * filename5 * ".jl" in list
    @test !(folder * "/" * filename4 * ".dat" in list)
    rm(base, recursive=true)
end