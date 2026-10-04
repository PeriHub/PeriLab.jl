# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const CF = PeriLab.InputDeck

cf_group(master, slave; extra...) = Dict{String,Any}("Master Block ID" => master,
                                                     "Slave Block ID" => slave,
                                                     "Search Radius" => 0.01,
                                                     (string(k) => v for (k, v) in extra)...)

function cf_contact(raw)
    ctx = PeriLab.ParameterSpec.ParseContext()
    return CF.parse_contact(raw, "Contact", ctx), ctx
end

function cf_section(T, raw, path = "test")
    ctx = PeriLab.ParameterSpec.ParseContext()
    return PeriLab.ParameterSpec.convert_value(T, raw, path, ctx), ctx
end

@testset "contact defaults" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("g" => cf_group(2, 1)))))
    @test isempty(ctx.errors)
    @test c.globals.global_search_frequency === 1
    @test c.globals.only_surface_contact_nodes
    m = c.models["C"]
    @test m.contact_stiffness === 1e8
    @test m.friction_coefficient === 0.0
    @test m.symmetry == "3D"
end

@testset "contact_search_frequency" begin
    c, ctx = cf_contact(Dict{String,Any}("Globals" => Dict{String,Any}("Global Search Frequency" => 3),
                                         "C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("own" => cf_group(2, 1; var"Global Search Frequency" = 5),
                                                                                                      "inherit" => cf_group(3, 4)))))
    @test isempty(ctx.errors)
    groups = c.models["C"].contact_groups
    @test CF.contact_search_frequency(groups["own"], c.globals) === 5
    @test CF.contact_search_frequency(groups["inherit"], c.globals) === 3
end

@testset "contact_blocks" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("a" => cf_group(1, 2),
                                                                                                      "c" => cf_group(3, 5),
                                                                                                      "q" => cf_group(8, 5))),
                                         "D" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("ba" => cf_group(2, 1)))))
    @test isempty(ctx.errors)
    @test CF.contact_blocks(c) == [1, 2, 3, 5, 8]
end

@testset "contact group checks" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("self" => cf_group(1, 1),
                                                                                                      "flat" => Dict{String,Any}("Master Block ID" => 1,
                                                                                                                                 "Slave Block ID" => 2,
                                                                                                                                 "Search Radius" => 0.0)))))
    @test length(ctx.errors) == 2
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["Contact.C.\"Contact Groups\".self"] ==
          "Master Block ID and Slave Block ID are equal; self contact is not implemented"
    @test msgs["Contact.C.\"Contact Groups\".flat.\"Search Radius\""] ==
          "must be greater than zero"
end

@testset "FEM degree" begin
    fem(degree) = Dict{String,Any}("Element Type" => "Lagrange", "Degree" => degree,
                                   "Material Model" => "Elastic")
    f, ctx = cf_section(CF.FEMParams, fem(2))
    @test isempty(ctx.errors) && CF.fem_degree(f) == [2]
    f, ctx = cf_section(CF.FEMParams, fem("2 1 1"))
    @test isempty(ctx.errors) && CF.fem_degree(f) == [2, 1, 1]
    for bad in ("2 a", "", "0 1", 0)
        f, ctx = cf_section(CF.FEMParams, fem(bad), "FEM")
        @test length(ctx.errors) == 1 && ctx.errors[1].path == "FEM.Degree"
    end
end

@testset "coupling defaults" begin
    f, ctx = cf_section(CF.FEMParams,
                        Dict{String,Any}("Element Type" => "Lagrange", "Degree" => 1,
                                         "Material Model" => "Elastic",
                                         "Coupling" => Dict{String,Any}("Coupling Type" => "Arlequin")))
    @test isempty(ctx.errors)
    @test f.coupling.pd_weight === 0.5
    @test f.coupling.kappa === 1.0
    @test f.coupling.coupling_block === nothing
end

@testset "surface correction type" begin
    s, ctx = cf_section(CF.SurfaceCorrectionParams,
                        Dict{String,Any}("Type" => "Area Correction"), "\"Surface Correction\"")
    @test length(ctx.errors) == 1 && ctx.errors[1].path == "\"Surface Correction\".Type"
end
