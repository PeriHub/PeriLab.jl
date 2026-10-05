# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const MPS = PeriLab.ParameterSpec
const MAT = PeriLab.Solver_Manager.Model_Factory.Material

function ut_material(dict; directory = "")
    ctx = MPS.ParseContext(directory = directory)
    m = MPS.parse_model(:material, Dict{String,Any}(dict), "Models.\"Material Models\".M",
                        ctx; name_key = "Material Model")
    return m, ctx
end

@testset "material base part" begin
    @test MPS.base_model(:material) === MAT.MaterialBaseParams
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic", "Symmetry" => "isotropic",
                              "Bulk Modulus" => 2.5e3, "Shear Modulus" => 1.15e3,
                              "C11" => 1.0, "State Factor ID" => 2,
                              "Zero Energy Control" => "Global", "Bond Associated" => true,
                              "Flaw Function" => Dict("Active" => true, "Function" => "Pre-defined",
                                                      "Flaw Size" => 0.2, "Flaw Magnitude" => 0.5,
                                                      "Flaw Location X" => 1.0,
                                                      "Flaw Location Y" => 0.0)))
    @test isempty(ctx.errors)
    @test m isa MPS.WithBase
    @test m.base.bulk_modulus == 2.5e3 && m.base.c11 == 1.0 && m.base.c66 === nothing
    @test m.base.bond_associated && !m.base.linear_strain
    @test m.base.flaw_function.flaw_size == 0.2
end

@testset "non-correspondence material names" begin
    for (name, T) in [("Bond-based Elastic", MAT.Bondbased_Elastic.BondbasedElasticParams),
                      ("1D Bond-based Elastic",
                       MAT.OneD_Bond_Based_Elastic.OneDBondbasedElasticParams),
                      ("Unified Bond-based Elastic",
                       MAT.Unified_Bondbased_Elastic.UnifiedBondbasedElasticParams),
                      ("PD Solid Elastic", MAT.PD_Solid_Elastic.PDSolidElasticParams),
                      ("PD Solid Plastic", MAT.PD_Solid_Plastic.PDSolidPlasticParams),
                      ("Rigid", MAT.Rigid.RigidParams)]
        @test MPS.lookup_model(:material, name) === T
    end
    m, ctx = ut_material(Dict("Material Model" => "1D Bond-based Elastic",
                              "Young's Modulus" => 1.0, "Id1" => 1, "Id2" => 2))
    @test isempty(ctx.errors) && m.model.id2 === 2
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic + PD Solid Plastic",
                              "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                              "Yield Stress" => 5))
    @test isempty(ctx.errors)
    @test m.model isa MPS.Composite
    @test m.model.parts[2].yield_stress == MPS.Constant(5.0)
end

@testset "PD Solid Plastic requires Yield Stress" begin
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Plastic", "Bulk Modulus" => 1.0))
    @test m === nothing
    @test only(ctx.errors).message == "missing (required by PD Solid Plastic)"
    @test only(ctx.errors).path == "Models.\"Material Models\".M.\"Yield Stress\""
end

@testset "material base constraints" begin
    m, ctx = ut_material(Dict("Material Model" => "Bond-based Elastic",
                              "Poisson's Ratio" => 0.7, "Shear Modulus" => -1.0,
                              "Flaw Function" => Dict("Active" => true, "Function" => "Gauss")))
    @test m === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    p = "Models.\"Material Models\".M"
    @test msgs["$p.\"Poisson's Ratio\""] == "0.7 is above maximum 0.5"
    @test msgs["$p.\"Shear Modulus\""] == "-1.0 is below minimum 0"
    @test msgs["$p.\"Flaw Function\".Function"] == "\"Gauss\" is not one of: \"Pre-defined\""
end

const CORR = MAT.Correspondence

@testset "correspondence material names" begin
    for (name, T) in [("Correspondence Elastic", CORR.Correspondence_Elastic.CorrespondenceElasticParams),
                      ("Correspondence Plastic", CORR.Correspondence_Plastic.CorrespondencePlasticParams),
                      ("Correspondence UMAT", CORR.Correspondence_UMAT.CorrespondenceUMATParams),
                      ("Correspondence VUMAT", CORR.Correspondence_VUMAT.CorrespondenceVUMATParams)]
        @test MPS.lookup_model(:material, name) === T
    end
    m, ctx = ut_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                              "Symmetry" => "isotropic plane strain", "Bulk Modulus" => 1.0,
                              "Shear Modulus" => 1.0, "Yield Stress" => 2.0))
    @test isempty(ctx.errors)
    @test m.model.parts[2].yield_stress == MPS.Constant(2.0)
end

@testset "UMAT properties" begin
    props = Dict("Property_$i" => Float64(i) for i in 1:27)
    m, ctx = ut_material(merge(Dict{String,Any}("Material Model" => "Correspondence UMAT",
                                                "File" => "libusertest.so",
                                                "Number of Properties" => 27,
                                                "Number of State Variables" => 0,
                                                "UMAT Material Name" => "test"), props))
    @test isempty(ctx.errors)
    @test m.model.number_of_properties == 27 && m.model.umat_name === nothing
    m, ctx = ut_material(Dict("Material Model" => "Correspondence VUMAT", "File" => "a.so",
                              "Number of Properties" => 3, "Property_3" => "abc",
                              "Propery_2" => 1.0))
    @test length(ctx.errors) == 2
    msgs = Dict(e.path => e.message for e in ctx.errors)
    p = "Models.\"Material Models\".M"
    @test msgs["$p.Property_3"] == "expected a number, got \"abc\""
    @test startswith(msgs["$p.Propery_2"], "unknown key")
    m, ctx = ut_material(Dict("Material Model" => "Correspondence UMAT", "File" => "a.so"))
    @test only(ctx.errors).message == "missing (required by Correspondence UMAT)"
end

@testset "the material template does not register (copies must not collide)" begin
    @test MPS.lookup_model(:material, "Material Template") === nothing
    @test MPS.parameter_spec(MAT.Material_template.MaterialTemplateParams) isa Vector
end

@testset "orthotropic completeness" begin
    full = Dict{String,Any}("Material Model" => "PD Solid Elastic", "Symmetry" => "Orthotropic",
                            "Young's Modulus X" => 1.0, "Young's Modulus Y" => 1.0,
                            "Young's Modulus Z" => 1.0, "Poisson's Ratio XY" => 0.3,
                            "Poisson's Ratio YZ" => 0.3, "Poisson's Ratio XZ" => 0.3,
                            "Shear Modulus XY" => 1.0, "Shear Modulus YZ" => 1.0,
                            "Shear Modulus XZ" => 1.0)
    m, ctx = ut_material(full)
    @test isempty(ctx.errors)
    delete!(full, "Shear Modulus XZ")
    m, ctx = ut_material(full)
    @test only(ctx.errors).path == "Models.\"Material Models\".M.Symmetry"
    @test only(ctx.errors).message == "\"Orthotropic\" requires Shear Modulus XZ"
end

@testset "anisotropic and transverse isotropic completeness" begin
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic",
                              "Symmetry" => "anisotropic", "C11" => 1.0))
    @test startswith(only(ctx.errors).message, "\"anisotropic\" requires C12, C13")
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic",
                              "Symmetry" => "transverse isotropic plane stress",
                              "Young's Modulus X" => 1.0, "Young's Modulus Y" => 1.0,
                              "Poisson's Ratio XY" => 0.3))
    @test only(ctx.errors).message ==
          "\"transverse isotropic plane stress\" requires Shear Modulus XY"
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic",
                              "Symmetry" => "isotropic plane strain",
                              "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    @test isempty(ctx.errors)
end

