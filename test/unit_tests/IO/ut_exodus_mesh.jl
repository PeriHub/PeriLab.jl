# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus

# Exodus mesh file with one block of `elem_type`; conn has one column per element
function exodus_mesh(coords::Matrix{Float64}, elem_type::String, conn::Matrix{Int32};
                     nodal = Dict{String,Vector{Float64}}())
    file = tempname() * ".e"
    init = Initialization{Int32(3),Int32(size(coords, 2)),Int32(size(conn, 2)),Int32(1),
                          Int32(0),Int32(0)}()
    exo = ExodusDatabase{Int32,Int32,Int32,Float64}(file, "w", init)
    write_coordinates(exo, coords)
    write_block(exo, 1, elem_type, conn)
    if !isempty(nodal)
        write_time(exo, 1, 0.0)
        write_names(exo, NodalVariable, collect(keys(nodal)))
        for (name, values) in nodal
            write_values(exo, NodalVariable, 1, name, values)
        end
    end
    close(exo)
    return file
end

function read_exodus_mesh(file)
    PeriLab.IO.Mesh_Import.Exodus_Mesh.read_mesh(typed_input(Dict()), file)
end

@testset "ut_exodus_mesh_tetra" begin
    # unit tetrahedron and a second one shifted by 2 in x
    coords = [0.0 1.0 0.0 0.0 2.0 3.0 2.0 2.0
              0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0
              0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0]
    file = exodus_mesh(coords, "TETRA", Int32[1 5; 2 6; 3 7; 4 8])
    mesh = read_exodus_mesh(file)
    rm(file)
    @test size(mesh, 1) == 2
    @test mesh.x ≈ [0.25, 2.25]
    @test mesh.y ≈ [0.25, 0.25]
    @test mesh.z ≈ [0.25, 0.25]
    @test mesh.volume ≈ [1 / 6, 1 / 6]
    @test mesh.block_id == [1, 1]
    # nodal variables of element meshes: mean over the element nodes
    file = exodus_mesh(coords, "TETRA", Int32[1 5; 2 6; 3 7; 4 8];
                       nodal = Dict("Temperature" =>
                                        [1.0, 2.0, 3.0, 4.0, 10.0, 10.0, 10.0,
                                            30.0]))
    mesh = read_exodus_mesh(file)
    rm(file)
    @test mesh.Temperature ≈ [2.5, 15.0]
end

@testset "ut_exodus_mesh_sphere" begin
    coords = [0.0 1.0 2.0
              0.0 0.5 0.0
              0.0 0.0 3.0]
    # nodal variables are taken over as mesh columns in element order (here 3, 1, 2)
    file = exodus_mesh(coords, "SPHERE", Int32[3 1 2];
                       nodal = Dict("Volume" => [0.1, 0.2, 0.3],
                                    "Temperature" => [300.0, 310.0, 320.0]))
    mesh = read_exodus_mesh(file)
    rm(file)
    @test mesh.x ≈ [2.0, 0.0, 1.0]
    @test mesh.z ≈ [3.0, 0.0, 0.0]
    @test mesh.volume ≈ [0.3, 0.1, 0.2]
    @test mesh.Volume ≈ [0.3, 0.1, 0.2]
    @test mesh.Temperature ≈ [320.0, 300.0, 310.0]
    @test mesh.block_id == [1, 1, 1]

    file = exodus_mesh(coords, "SPHERE", Int32[1 2 3];
                       nodal = Dict("Temperature" => [300.0, 310.0, 320.0]))
    @test_logs (:error,
                "Volume is missing. Please define a 'Volume' for each point in the mesh file.") match_mode=:any @test_throws PeriLab.PeriLabError read_exodus_mesh(file)
    rm(file)
end

@testset "ut_exodus_mesh_unsupported" begin
    file = exodus_mesh([0.0 1.0; 0.0 0.0; 0.0 0.0], "BAR2", reshape(Int32[1, 2], 2, 1))
    @test_logs (:error, "Element type BAR2 not supported") match_mode=:any @test_throws PeriLab.PeriLabError read_exodus_mesh(file)
    rm(file)
end
