# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const BMAT = PeriLab.Solver_Manager.Model_Factory.Material
const BBASIS = PeriLab.Solver_Manager.Material_Basis

function ut_reset(dof; nnodes = 3)
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(nnodes)
    PeriLab.Data_Manager.set_dof(dof)
end

# legacy completion and typed completion on the same raw block, each on a fresh Data_Manager
function ut_both(raw; dof = 3, setup = () -> nothing)
    ut_reset(dof)
    setup()
    legacy = Dict{String,Any}(raw)
    BBASIS.get_all_elastic_moduli(legacy)
    legacy_values = Dict(k => copy(legacy[k])
                         for k in ("Bulk Modulus", "Young's Modulus", "Shear Modulus",
                                   "Poisson's Ratio"))
    ut_reset(dof)
    setup()
    typed = typed_block_material(raw; dof = dof)
    return legacy_values, typed
end

function ut_same(legacy, m::BMAT.ElasticModuli)
    return isapprox(legacy["Bulk Modulus"], m.bulk_modulus) &&
           isapprox(legacy["Young's Modulus"], m.youngs_modulus) &&
           isapprox(legacy["Shear Modulus"], m.shear_modulus) &&
           isapprox(legacy["Poisson's Ratio"], m.poissons_ratio)
end

@testset "isotropic completion matches the legacy completion" begin
    for raw in [Dict("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 10.0,
                     "Shear Modulus" => 10.0),
                Dict("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 5.0,
                     "Young's Modulus" => 1.25),
                Dict("Material Model" => "PD Solid Elastic", "Poisson's Ratio" => 0.45,
                     "Shear Modulus" => 1.25),
                Dict("Material Model" => "PD Solid Elastic", "Young's Modulus" => 5.0,
                     "Poisson's Ratio" => 0.125, "Symmetry" => "isotropic"),
                Dict("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 1.0,
                     "Shear Modulus" => 10.0, "Poisson's Ratio" => 0.2),
                Dict("Material Model" => "Unified Bond-based Elastic",
                     "Young's Modulus" => 5.0, "Poisson's Ratio" => 0.125,
                     "Symmetry" => "isotropic plane strain")]
        legacy, typed = ut_both(raw)
        @test typed isa BMAT.BlockMaterial
        @test ut_same(legacy, typed.moduli)
    end
end

@testset "bond-based fixed Poisson's ratio" begin
    for (dof, nu) in ((2, 1 / 3), (3, 1 / 4))
        legacy, typed = ut_both(Dict("Material Model" => "Bond-based Elastic",
                                     "Young's Modulus" => 5.0, "Poisson's Ratio" => 0.125,
                                     "Symmetry" => dof == 2 ? "isotropic plane stress" :
                                                   "isotropic"); dof = dof)
        @test typed.moduli.poissons_ratio == nu
        @test typed.moduli.youngs_modulus == 5.0
        @test ut_same(legacy, typed.moduli)
    end
end

@testset "moduli from a mesh field" begin
    setup = () -> PeriLab.Data_Manager.create_constant_node_scalar_field("Bulk_Modulus",
                                                                        Float64;
                                                                        default_value = 10)
    legacy, typed = ut_both(Dict("Material Model" => "PD Solid Elastic",
                                 "Shear Modulus" => 10.0); setup = setup)
    @test typed.moduli.youngs_modulus == [22.5, 22.5, 22.5]
    @test ut_same(legacy, typed.moduli)
    @test PeriLab.Data_Manager.get_field("Young's_Modulus") == [22.5, 22.5, 22.5]
    @test BMAT.modulus(typed.moduli.youngs_modulus, 2) == 22.5
    @test BMAT.modulus(5.0, 2) == 5.0
end

@testset "Hooke-matrix symmetries have no isotropic moduli" begin
    ut_reset(3)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Symmetry" => "orthotropic",
                                  "Young's Modulus X" => 1.0, "Young's Modulus Y" => 1.0,
                                  "Young's Modulus Z" => 1.0, "Poisson's Ratio XY" => 0.3,
                                  "Poisson's Ratio YZ" => 0.3, "Poisson's Ratio XZ" => 0.3,
                                  "Shear Modulus XY" => 1.0, "Shear Modulus YZ" => 1.0,
                                  "Shear Modulus XZ" => 1.0))
    @test m.moduli === nothing
end

@testset "too few isotropic constants" begin
    ut_reset(3)
    @test_logs (:error,
                "Minimum of two parameters are needed for isotropic material") @test_throws PeriLab.PeriLabError begin
        typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 10.0))
    end
end

@testset "material_symmetry follows check_symmetry and get_symmetry" begin
    for (sym, dof) in (("isotropic plane strain", 2), ("isotropic plane stress", 2),
                       ("isotropic plane strain", 3), ("Isotropic Plane Stress", 2),
                       ("isotropic", 3), (nothing, 2), (nothing, 3),
                       ("iso Plane Strain", 3))
        legacy = sym === nothing ? Dict{String,Any}() : Dict{String,Any}("Symmetry" => sym)
        if dof == 3 && haskey(legacy, "Symmetry")
            legacy["Symmetry"] = replace(replace(legacy["Symmetry"], r"plane strain$" => ""),
                                         r"plane stress$" => "")
        end
        @test BMAT.material_symmetry(sym, dof) == BBASIS.get_symmetry(legacy)
    end
end

@testset "write_moduli! fills the material dict" begin
    ut_reset(3)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 10.0, "Shear Modulus" => 10.0))
    dict = Dict{String,Any}("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 10.0,
                            "Shear Modulus" => 10.0)
    BMAT.write_moduli!(dict, m)
    @test dict["Young's Modulus"] == 22.5 && dict["Poisson's Ratio"] == 0.125
    @test dict["Computed"] === true
    @test dict["Symmetry"] == "isotropic"
    @test BMAT.model_parts(m.model) == (m.model,)
end

@testset "factory dispatches typed parts to their modules" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "Bond-based Elastic",
                                  "Young's Modulus" => 1.0))
    @test parentmodule(typeof(m.model)) === BMAT.Bondbased_Elastic
    @test hasmethod(BMAT.Bondbased_Elastic.compute_model,
                    Tuple{Vector{Int64},typeof(m.model),typeof(m),Int64,Float64,Float64})
    @test hasmethod(BMAT.Rigid.init_model,
                    Tuple{Vector{Int64},BMAT.Rigid.RigidParams,Any})
end

@testset "table follows the NP1 switch" begin
    ut_reset(3; nnodes = 2)
    dir = mktempdir()
    write(joinpath(dir, "ys.txt"), "header: Temperature Yield_Stress\n0 10\n100 20\n")
    N, NP1 = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    NP1 .= [0.0, 100.0]
    ctx = PeriLab.ParameterSpec.ParseContext(directory = dir)
    wb = PeriLab.ParameterSpec.parse_model(:material,
                                           Dict{String,Any}("Material Model" => "PD Solid Plastic",
                                                            "Bulk Modulus" => 1.0,
                                                            "Shear Modulus" => 1.0,
                                                            "Yield Stress" => "ys.txt"),
                                           "M", ctx; name_key = "Material Model")
    @test isempty(ctx.errors)
    m = BMAT.block_material(wb, "PD Solid Plastic", 3)
    BMAT.bind_material!(m)
    @test PeriLab.ParameterSpec.value(m.model.yield_stress, 2) ≈ 20.0
    PeriLab.Data_Manager.create_constant_node_scalar_field("Active", Bool; default_value = true)
    PeriLab.Data_Manager.switch_NP1_to_N()
    PeriLab.Data_Manager.get_field("Temperature", "NP1") .= [100.0, 0.0]
    BMAT.bind_material!(m)
    @test PeriLab.ParameterSpec.value(m.model.yield_stress, 1) ≈ 20.0
    @test PeriLab.ParameterSpec.value(m.model.yield_stress, 2) ≈ 10.0
end

@testset "composite parts in order" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic + PD Solid Plastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                  "Yield Stress" => 2.0))
    @test parentmodule.(typeof.(BMAT.model_parts(m.model))) ==
          (BMAT.PD_Solid_Elastic, BMAT.PD_Solid_Plastic)
end

