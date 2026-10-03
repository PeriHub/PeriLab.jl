# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_section(T, dict)
    ctx = PS.ParseContext()
    return PS.parse_section(T, dict, "X", ctx), ctx
end

@testset "Discretization: minimal and full" begin
    d, ctx = ut_section(ID.DiscretizationParams,
                        Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "mesh.txt"))
    @test isempty(ctx.errors)
    @test d.type == "Text File" && d.input_mesh_file == "mesh.txt"
    @test isempty(d.node_sets) && isempty(d.bond_filters)
    @test d.gcode === nothing && d.surface_extrusion === nothing
    @test d.influence_function === nothing
    full = Dict{String,Any}("Type" => "Exodus", "Input Mesh File" => "m.g",
                            "Distribution Type" => "Neighbor based",
                            "Influence Function" => "1/xi^2",
                            "Horizon Mesh Scaling X" => 1.5,
                            "Input External Topology" => Dict{String,Any}("File" => "t.txt",
                                                                          "Add Neighbor Search" => true),
                            "Surface Extrusion" => Dict{String,Any}("Direction" => "X",
                                                                    "Step_X" => 0.1,
                                                                    "Step_Y" => 0.1,
                                                                    "Step_Z" => 0,
                                                                    "Number" => 3),
                            "Gcode" => Dict{String,Any}("Overwrite Mesh" => true,
                                                        "Sampling" => 0.5, "Width" => 1,
                                                        "Height" => 0.2,
                                                        "Blocks" => Dict{String,Any}("1" => "block_1")),
                            "Bond Filters" => Dict{String,Any}("bf_1" => Dict{String,Any}("Type" => "Rectangular_Plane",
                                                                                          "Normal X" => 0.0,
                                                                                          "Normal Y" => 1.0,
                                                                                          "Allow Contact" => true)))
    d, ctx = ut_section(ID.DiscretizationParams, full)
    @test isempty(ctx.errors)
    @test d.influence_function == "1/xi^2"
    @test d.horizon_mesh_scaling_x === 1.5 && d.horizon_mesh_scaling_y === nothing
    @test d.input_external_topology.add_neighbor_search === true
    @test d.surface_extrusion.step_z === 0.0 && d.surface_extrusion.number === 3
    @test d.gcode.scale === 1.0 && d.gcode.blocks == Dict("1" => "block_1")
    @test d.bond_filters["bf_1"].allow_contact && d.bond_filters["bf_1"].normal_z === nothing
end

@testset "node set entries" begin
    d, ctx = ut_section(ID.DiscretizationParams,
                        Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                                         "Node Sets" => Dict{String,Any}("a" => 5,
                                                                         "b" => "2 3 4",
                                                                         "c" => "ns10.txt")))
    @test isempty(ctx.errors)
    @test d.node_sets["a"] === 5
    @test d.node_sets["b"] == "2 3 4"
    @test d.node_sets["c"] == "ns10.txt"
end

@testset "Discretization errors" begin
    d, ctx = ut_section(ID.DiscretizationParams,
                        Dict{String,Any}("Type" => "Text File",
                                         "Surface Extrusion" => Dict{String,Any}("Direction" => "W",
                                                                                 "Step_X" => 1,
                                                                                 "Step_Y" => 1,
                                                                                 "Step_Z" => 1,
                                                                                 "Number" => 1)))
    @test d === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["X.\"Input Mesh File\""] == "missing (required by DiscretizationParams)"
    @test msgs["X.\"Surface Extrusion\".Direction"] ==
          "\"W\" is not one of: \"X\", \"Y\", \"Z\""
end

@testset "Blocks" begin
    b, ctx = ut_section(ID.BlockParams,
                        Dict{String,Any}("Block ID" => 1, "Density" => 2700,
                                         "Horizon" => 2, "Material Model" => "Steel",
                                         "Degradation Model" => "deg", "FEM" => true,
                                         "Step ID" => "1,2"))
    @test isempty(ctx.errors)
    @test b.block_id === 1 && b.density === 2700.0 && b.horizon === 2.0
    @test b.material_model == "Steel" && b.damage_model === nothing
    @test b.degradation_model == "deg" && b.fem === true && b.step_id == "1,2"
    b, ctx = ut_section(ID.BlockParams,
                        Dict{String,Any}("Block ID" => 1, "Density" => 1.0, "Horizon" => -1.0))
    @test ctx.errors[1].message == "-1.0 is below minimum 0"
end

@testset "FEM and Surface Correction" begin
    f, ctx = ut_section(ID.FEMParams,
                        Dict{String,Any}("Element Type" => "Lagrange", "Degree" => "1 1",
                                         "Material Model" => "Elastic",
                                         "Coupling" => Dict{String,Any}("Coupling Type" => "Arlequin",
                                                                        "PD Weight" => 0.5,
                                                                        "Coupling Block" => 2)))
    @test isempty(ctx.errors)
    @test f.degree == "1 1" && f.coupling.pd_weight === 0.5 && f.coupling.coupling_block === 2
    f, ctx = ut_section(ID.FEMParams,
                        Dict{String,Any}("Element Type" => "Lagrange", "Degree" => 1,
                                         "Material Model" => "Elastic"))
    @test f.degree === 1 && f.coupling === nothing
    s, ctx = ut_section(ID.SurfaceCorrectionParams,
                        Dict{String,Any}("Type" => "Volume Correction"))
    @test isempty(ctx.errors) && s.update === false
end
