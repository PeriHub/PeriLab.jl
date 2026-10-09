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
    @test m.name == "Penalty Contact"
    @test m.model.contact_stiffness === 1e8
    @test m.model.friction_coefficient === 0.0
    @test m.base.symmetry == "3D"
    @test m.base.contact_radius === 0.005
end

# a contact model declares its own keys; the shared keys come from the base part
PeriLab.ParameterSpec.@params struct UTContactModelParams
    ut_value::Float64 = req("UT Value"; min = 0)
end
PeriLab.ParameterSpec.register_contact("UT Contact", UTContactModelParams)

cf_model(type; extra...) = Dict{String,Any}("Type" => type, "Contact Radius" => 0.005,
                                            "Contact Groups" => Dict{String,Any}("g" => cf_group(2, 1)),
                                            (string(k) => v for (k, v) in extra)...)

@testset "contact models declare their own keys" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => cf_model("UT Contact"; var"UT Value" = 2.0)))
    @test isempty(ctx.errors)
    m = c.models["C"]
    @test m.model isa UTContactModelParams && m.model.ut_value === 2.0
    @test m.base.contact_groups["g"].master_block_id === 2
    @test CF.contact_blocks(c) == [1, 2]
    # Penalty keys belong to the penalty model only
    c, ctx = cf_contact(Dict{String,Any}("C" => cf_model("UT Contact"; var"UT Value" = 2.0,
                                                         var"Contact Stiffness" = 1.0)))
    @test ctx.errors[1].path == "Contact.C.\"Contact Stiffness\""
    @test startswith(ctx.errors[1].message, "unknown key")
end

@testset "contact model names" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => cf_model("Penalty Contakt")))
    @test c === nothing
    @test ctx.errors[1].path == "Contact.C.Type"
    @test ctx.errors[1].message ==
          "model \"Penalty Contakt\" not found — did you mean \"Penalty Contact\"?"
    c, ctx = cf_contact(Dict{String,Any}("C" => cf_model("Penalty Contact + UT Contact";
                                                         var"UT Value" = 2.0)))
    @test c === nothing
    @test ctx.errors[1].path == "Contact.C.Type"
    @test ctx.errors[1].message == "contact models cannot be combined with +"
end

@testset "contact_search_frequency" begin
    c, ctx = cf_contact(Dict{String,Any}("Globals" => Dict{String,Any}("Global Search Frequency" => 3),
                                         "C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("own" => cf_group(2, 1; var"Global Search Frequency" = 5),
                                                                                                      "inherit" => cf_group(3, 4)))))
    @test isempty(ctx.errors)
    groups = c.models["C"].base.contact_groups
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
