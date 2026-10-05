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


@testset "binding is cheap when a material has no tables" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    BMAT.bind_material!(m)
    @test (@allocated BMAT.bind_material!(m)) < 1000
end

@testset "block materials are not written to checkpoints" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    PeriLab.Data_Manager.set_block_material(1, m)
    dir = mktempdir()
    PeriLab.Data_Manager.write_checkpoint(dir, 0)
    PeriLab.Data_Manager.set_block_material(1, nothing)
    PeriLab.Data_Manager.read_checkpoint!(dir, 0)
    # rebuilt from the input by read_properties, never restored from a checkpoint
    @test PeriLab.Data_Manager.get_block_material(1) === nothing
end

function ut_legacy_dict(raw; dof)
    ut_reset(dof)
    legacy = Dict{String,Any}(raw)
    BBASIS.get_all_elastic_moduli(legacy)
    return legacy
end

@testset "typed Hooke matrix equals legacy" begin
    ortho = Dict("Young's Modulus X" => 2.0, "Young's Modulus Y" => 1.5,
                 "Young's Modulus Z" => 1.0, "Poisson's Ratio XY" => 0.3,
                 "Poisson's Ratio YZ" => 0.25, "Poisson's Ratio XZ" => 0.2,
                 "Shear Modulus XY" => 0.7, "Shear Modulus YZ" => 0.6,
                 "Shear Modulus XZ" => 0.5)
    aniso = Dict("C$i$j" => (i == j ? 100.0 * i : 1.0 * i + j) for i in 1:6 for j in i:6)
    cases = [(3, Dict("Symmetry" => "isotropic", "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0)),
             (2, Dict("Symmetry" => "isotropic plane strain", "Bulk Modulus" => 10.0,
                      "Shear Modulus" => 4.0)),
             (2, Dict("Symmetry" => "isotropic plane stress", "Young's Modulus" => 10.0,
                      "Poisson's Ratio" => 0.3)),
             (2, Dict("Symmetry" => "something else", "Young's Modulus" => 10.0,
                      "Poisson's Ratio" => 0.3)),
             (3, merge(Dict("Symmetry" => "orthotropic"), ortho)),
             (2, merge(Dict("Symmetry" => "orthotropic plane strain"), ortho)),
             (3, merge(Dict("Symmetry" => "transverse isotropic"), ortho)),
             (2, merge(Dict("Symmetry" => "transverse isotropic plane strain"), ortho)),
             (2, merge(Dict("Symmetry" => "transverse isotropic plane stress"), ortho)),
             (3, merge(Dict{String,Any}("Symmetry" => "anisotropic"), aniso)),
             (2, merge(Dict{String,Any}("Symmetry" => "anisotropic plane strain"), aniso))]
    for (dof, raw) in cases
        raw = merge(Dict{String,Any}("Material Model" => "Correspondence Elastic"), raw)
        legacy = ut_legacy_dict(raw; dof = dof)
        ut_reset(dof)
        m = typed_block_material(raw; dof = dof)
        @test Matrix(BBASIS.hooke_matrix(m, dof, 2)) ≈
              Matrix(BBASIS.get_Hooke_matrix(legacy, legacy["Symmetry"], dof, 2))
    end
end

@testset "hooke_symmetry follows the dict the legacy path used" begin
    @test BMAT.hooke_symmetry(nothing, 3) == "isotropic"
    @test BMAT.hooke_symmetry("isotropic plane strain", 3) == "isotropic "
    @test BMAT.hooke_symmetry("isotropic plane strain", 2) == "isotropic plane strain"
    @test BMAT.hooke_symmetry("orthotropic", 3) == "orthotropic"
end

@testset "Hooke matrix from a table" begin
    ut_reset(3; nnodes = 2)
    dir = mktempdir()
    write(joinpath(dir, "ex.txt"), "header: Temperature Young's_Modulus_X\n0 1000\n100 3000\n")
    N, NP1 = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    NP1 .= [0.0, 100.0]
    raw = Dict{String,Any}("Material Model" => "Correspondence Elastic",
                           "Symmetry" => "orthotropic", "Young's Modulus X" => "ex.txt",
                           "Young's Modulus Y" => 1.5e3, "Young's Modulus Z" => 1.0e3,
                           "Poisson's Ratio XY" => 0.3, "Poisson's Ratio YZ" => 0.25,
                           "Poisson's Ratio XZ" => 0.2, "Shear Modulus XY" => 700.0,
                           "Shear Modulus YZ" => 600.0, "Shear Modulus XZ" => 500.0)
    ctx = PeriLab.ParameterSpec.ParseContext(directory = dir)
    wb = PeriLab.ParameterSpec.parse_model(:material, raw, "M", ctx;
                                           name_key = "Material Model")
    @test isempty(ctx.errors)
    m = BMAT.block_material(wb, "Correspondence Elastic", 3)
    BMAT.bind_material!(m)
    constant_x(E) = begin
        r = copy(raw)
        r["Young's Modulus X"] = E
        ut_reset(3; nnodes = 2)
        BBASIS.hooke_matrix(typed_block_material(r; dof = 3), 3, 1)
    end
    @test Matrix(BBASIS.hooke_matrix(m, 3, 1)) ≈ Matrix(constant_x(1000.0))
    @test Matrix(BBASIS.hooke_matrix(m, 3, 2)) ≈ Matrix(constant_x(3000.0))
end

@testset "typed flaw function" begin
    ut_reset(3)
    flaw = Dict("Active" => true, "Function" => "Pre-defined", "Flaw Size" => 0.2,
                "Flaw Magnitude" => 0.5, "Flaw Location X" => 1.0, "Flaw Location Y" => 0.5)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                  "Flaw Function" => flaw))
    for coor in ([1.0, 0.5], [1.1, 0.4, 0.2], [3.0, 3.0])
        @test BBASIS.flaw_function(m.base.flaw_function, coor, 10.0) ≈
              BBASIS.flaw_function(Dict("Flaw Function" => flaw), coor, 10.0)
    end
    @test BBASIS.flaw_function(nothing, [0.0, 0.0], 10.0) == 10.0
    inactive = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                         "Flaw Function" => Dict("Active" => false,
                                                                 "Function" => "Pre-defined")))
    @test BBASIS.flaw_function(inactive.base.flaw_function, [0.0, 0.0], 10.0) == 10.0
    incomplete = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                           "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                           "Flaw Function" => Dict("Active" => true,
                                                                   "Function" => "Pre-defined")))
    @test_logs (:error,
                "An active Flaw Function needs Flaw Size, Flaw Magnitude, Flaw Location X and Flaw Location Y.") @test_throws PeriLab.PeriLabError begin
        BBASIS.flaw_function(incomplete.base.flaw_function, [0.0, 0.0], 10.0)
    end
end

@testset "extras reach the block material" begin
    ut_reset(3)
    m = typed_block_material(Dict("Material Model" => "Correspondence UMAT", "File" => "x.so",
                                  "Number of Properties" => 2, "Property_1" => 3.0,
                                  "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    @test m.extras == Dict{String,Any}("Property_1" => 3.0)
    @test m.hooke_symmetry == "isotropic"
end


const BCORR = BMAT.Correspondence

@testset "Correspondence Elastic typed init equals legacy" begin
    for (dof, raw) in ((3, Dict("Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                "Shear Modulus" => 4.0)),
                       (2, Dict("Symmetry" => "isotropic plane strain",
                                "Young's Modulus" => 10.0, "Poisson's Ratio" => 0.3)))
        raw = merge(Dict{String,Any}("Material Model" => "Correspondence Elastic"), raw)
        legacy = ut_legacy_dict(raw; dof = dof)
        ut_reset(dof; nnodes = 2)
        m = typed_block_material(raw; dof = dof)
        BCORR.Correspondence_Elastic.init_model([1, 2], m.model, m)
        C = PeriLab.Data_Manager.get_field("Material Gradient")
        for iID in 1:2
            @test C[iID, :, :] ≈ Matrix(BBASIS.get_Hooke_matrix(legacy, legacy["Symmetry"],
                                                                dof, iID))
        end
        @test hasmethod(BCORR.Correspondence_Elastic.compute_stresses,
                        Tuple{Vector{Int64},Int64,typeof(m.model),typeof(m),Float64,Float64,
                              Array{Float64,3},Array{Float64,3},Array{Float64,3}})
    end
end

@testset "Correspondence Plastic typed compute equals legacy" begin
    dof = 3
    nnodes = 2
    raw = Dict{String,Any}("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                           "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                           "Shear Modulus" => 4.0, "Yield Stress" => 0.01)
    results = []
    for typed in (false, true)
        ut_reset(dof; nnodes = nnodes)
        coor = PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates", Float64,
                                                                      dof)
        coor .= [0.0 0.0 0.0; 1.0 0.0 0.0]
        if typed
            m = typed_block_material(raw; dof = dof)
            p = m.model.parts[2]
            BCORR.Correspondence_Plastic.init_model(collect(1:nnodes), p, m)
        else
            legacy = Dict{String,Any}(raw)
            BBASIS.get_all_elastic_moduli(legacy)
            BCORR.Correspondence_Plastic.init_model(collect(1:nnodes), legacy)
        end
        strain_inc = zeros(nnodes, dof, dof)
        strain_inc[:, 1, 1] .= 0.02
        strain_inc[:, 1, 2] .= 0.01
        strain_inc[:, 2, 1] .= 0.01
        stress_N = zeros(nnodes, dof, dof)
        stress_NP1 = zeros(nnodes, dof, dof)
        stress_NP1[:, 1, 1] .= 0.5
        if typed
            BCORR.Correspondence_Plastic.compute_stresses(collect(1:nnodes), dof, p, m, 0.0,
                                                          1.0, strain_inc, stress_N,
                                                          stress_NP1)
        else
            BCORR.Correspondence_Plastic.compute_stresses(collect(1:nnodes), dof, legacy,
                                                          0.0, 1.0, strain_inc, stress_N,
                                                          stress_NP1)
        end
        push!(results, copy(stress_NP1))
    end
    @test results[1] ≈ results[2]
end

@testset "zero energy control skips UMAT" begin
    ut_reset(3; nnodes = 2)
    elastic = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                        "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    umat = typed_block_material(Dict("Material Model" => "Correspondence UMAT",
                                     "File" => "x.so", "Number of Properties" => 1,
                                     "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    vumat = typed_block_material(Dict("Material Model" => "Correspondence VUMAT",
                                      "File" => "x.so", "Number of Properties" => 1,
                                      "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    GZEC = PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control
    @test !GZEC.is_umat(elastic)
    @test GZEC.is_umat(umat)
    @test !GZEC.is_umat(vumat)
end

@testset "typed zero energy control init" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0,
                                  "Zero Energy Control" => "Global"))
    ZEC = PeriLab.Solver_Manager.Zero_Energy_Control
    ZEC.init_model([1, 2], m, 1)
    @test PeriLab.Data_Manager.get_analysis_model("Zero Energy Control Model", 1) == ["Global"]
    @test PeriLab.Data_Manager.get_field("Material Gradient")[2, :, :] ≈
          Matrix(BBASIS.hooke_matrix(m, 3, 2))
    ut_reset(3; nnodes = 2)
    plain = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                      "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0))
    ZEC.init_model([1, 2], plain, 1)
    @test PeriLab.Data_Manager.get_analysis_model("Zero Energy Control Model", 1) == [""]
end


function ut_umat_file()
    file = "./src/Models/Material/UMATs/libperuser.so"
    isfile(file) || (file = "../src/Models/Material/UMATs/libperuser.so")
    return file
end

@testset "UMAT properties from extras" begin
    ut_reset(3; nnodes = 2)
    UMAT = BCORR.Correspondence_UMAT
    file = ut_umat_file()
    m = typed_block_material(Dict("Material Model" => "Correspondence UMAT", "File" => file,
                                  "Number of Properties" => 3, "Property_1" => 2,
                                  "Property_3" => 2.4, "Young's Modulus" => 2.0,
                                  "Poisson's Ratio" => 0.1))
    UMAT.init_model([1, 2], m.model, m)
    props = PeriLab.Data_Manager.get_field("Properties")
    @test props[1] == 2.0 && props[2] == 0.0 && props[3] == 2.4
    @test PeriLab.Data_Manager.get_field("Material Gradient")[1, :, :] ≈
          Matrix(BBASIS.hooke_matrix(m, 3, 1))
    @test UMAT.umat_file_path == joinpath(pwd(), PeriLab.Data_Manager.get_directory(), file)
end

@testset "UMAT and VUMAT typed init errors" begin
    for (UM, name_key, model) in ((BCORR.Correspondence_UMAT, "UMAT Material Name",
                                   "Correspondence UMAT"),
                                  (BCORR.Correspondence_VUMAT, "VUMAT Material Name",
                                   "Correspondence VUMAT"))
        ut_reset(3; nnodes = 2)
        file = ut_umat_file()
        missing_file = typed_block_material(Dict("Material Model" => model,
                                                 "File" => file * "_not_there",
                                                 "Number of Properties" => 1,
                                                 "Young's Modulus" => 2.0,
                                                 "Poisson's Ratio" => 0.1))
        @test_logs (:error,
                    "File $(joinpath(pwd(), PeriLab.Data_Manager.get_directory(), file * "_not_there")) does not exist, please check name and directory.") @test_throws PeriLab.PeriLabError begin
            UM.init_model([1, 2], missing_file.model, missing_file)
        end
        long_name = typed_block_material(Dict("Material Model" => model, "File" => file,
                                              "Number of Properties" => 1,
                                              name_key => "a"^81,
                                              "Young's Modulus" => 2.0,
                                              "Poisson's Ratio" => 0.1))
        @test_logs (:error,
                    "Due to old Fortran standards only a name length of 80 is supported") @test_throws PeriLab.PeriLabError begin
            UM.init_model([1, 2], long_name.model, long_name)
        end
    end
end

@testset "UMAT predefined fields (typed)" begin
    ut_reset(3; nnodes = 2)
    t2 = PeriLab.Data_Manager.create_constant_node_scalar_field("test_field_2", Float64)
    t2[1] = 7.3
    t3 = PeriLab.Data_Manager.create_constant_node_scalar_field("test_field_3", Float64)
    t3 .= 3
    m = typed_block_material(Dict("Material Model" => "Correspondence UMAT",
                                  "File" => ut_umat_file(), "Number of Properties" => 1,
                                  "Predefined Field Names" => "test_field_2 test_field_3",
                                  "Young's Modulus" => 2.0, "Poisson's Ratio" => 0.1))
    BCORR.Correspondence_UMAT.init_model([1, 2], m.model, m)
    fields = PeriLab.Data_Manager.get_field("Predefined Fields")
    @test fields[1, 1] == 7.3 && fields[2, 2] == 3.0
end


@testset "typed correspondence dispatcher" begin
    ut_reset(3; nnodes = 2)
    PeriLab.Data_Manager.set_rotation(false)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                  "Shear Modulus" => 4.0))
    BCORR.init_model([1, 2], 1, m)
    @test PeriLab.Data_Manager.get_analysis_model("Correspondence Model", 1) ==
          ["Correspondence Elastic"]
    @test PeriLab.Data_Manager.get_field("Material Gradient")[1, :, :] ≈
          Matrix(BBASIS.hooke_matrix(m, 3, 1))
    ut_reset(3; nnodes = 2)
    nosym = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                      "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0))
    # as before (the dict got "isotropic" when Symmetry was missing): 3D runs isotropic
    BCORR.init_model([1, 2], 1, nosym)
    @test PeriLab.Data_Manager.get_field("Material Gradient")[1, :, :] ≈
          Matrix(BBASIS.hooke_matrix(nosym, 3, 1))
    @test nosym.hooke_symmetry == "isotropic"
    @test hasmethod(BCORR.compute_model, Tuple{Vector{Int64},typeof(m),Int64,Float64,Float64})
    @test hasmethod(BCORR.Bond_Associated_Correspondence.compute_model,
                    Tuple{Vector{Int64},typeof(m),Int64,Float64,Float64})
end


@testset "UMAT state factor scales the moduli once" begin
    ut_reset(3; nnodes = 3)
    sv = PeriLab.Data_Manager.create_constant_node_scalar_field("State Variables", Float64)
    sv .= 2.0
    m = typed_block_material(Dict("Material Model" => "Correspondence UMAT",
                                  "File" => ut_umat_file(), "Number of Properties" => 1,
                                  "State Factor ID" => 1, "Young's Modulus" => 10.0,
                                  "Poisson's Ratio" => 0.25))
    E = copy(m.moduli.youngs_modulus)
    BCORR.Correspondence_UMAT.init_model([1, 2, 3], m.model, m)
    @test PeriLab.Data_Manager.get_field("Young's_Modulus") ≈ 2 .* E
end

@testset "UMAT init with a state factor grows linearly with the node count" begin
    function umat_init_allocations(nnodes)
        ut_reset(3; nnodes = nnodes)
        m = typed_block_material(Dict("Material Model" => "Correspondence UMAT",
                                      "File" => ut_umat_file(), "Number of Properties" => 1,
                                      "State Factor ID" => 1, "Young's Modulus" => 10.0,
                                      "Poisson's Ratio" => 0.25))
        nodes = collect(1:nnodes)
        return @allocated BCORR.Correspondence_UMAT.init_model(nodes, m.model, m)
    end
    umat_init_allocations(50)
    a1 = umat_init_allocations(500)
    a2 = umat_init_allocations(1000)
    @test a2 < 3 * a1
end

@testset "correspondence flag" begin
    ut_reset(3)
    @test typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                    "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0)).correspondence
    @test !typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                     "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0)).correspondence
end

@testset "critical bulk modulus" begin
    ut_reset(3)
    iso = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                    "Bulk Modulus" => 7.0, "Shear Modulus" => 2.0))
    @test BMAT.critical_bulk_modulus(iso) == 7.0
    ortho = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                      "Symmetry" => "orthotropic",
                                      "Young's Modulus X" => 2.0, "Young's Modulus Y" => 1.5,
                                      "Young's Modulus Z" => 1.0, "Poisson's Ratio XY" => 0.3,
                                      "Poisson's Ratio YZ" => 0.25, "Poisson's Ratio XZ" => 0.2,
                                      "Shear Modulus XY" => 0.7, "Shear Modulus YZ" => 0.6,
                                      "Shear Modulus XZ" => 0.5))
    s11, s22, s33 = 1 / 2.0, 1 / 1.5, 1 / 1.0
    s12, s23, s13 = -0.3 / 2.0, -0.25 / 1.0, -0.2 / 1.0
    @test BMAT.critical_bulk_modulus(ortho) ≈
          1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s13))
    aniso = typed_block_material(merge(Dict{String,Any}("Material Model" => "Correspondence Elastic",
                                                        "Symmetry" => "anisotropic"),
                                       Dict("C$i$j" => (i == j ? 10.0 * i : 1.0)
                                            for i in 1:6 for j in i:6)))
    @test BMAT.critical_bulk_modulus(aniso) == maximum([40.0 / 2, 50.0 / 2, 60.0 / 2])
    transverse = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                           "Symmetry" => "transverse isotropic plane stress",
                                           "Young's Modulus X" => 2.0,
                                           "Young's Modulus Y" => 1.5,
                                           "Poisson's Ratio XY" => 0.3,
                                           "Shear Modulus XY" => 0.8))
    @test BMAT.critical_bulk_modulus(transverse) == 0.4
end

@testset "pre-calculation dependencies" begin
    function ut_dependencies(raw)
        ut_reset(3; nnodes = 2)
        PeriLab.Data_Manager.set_block_id_list([1])
        PeriLab.Data_Manager.init_properties()
        PeriLab.Data_Manager.set_block_material(1, typed_block_material(raw))
        PeriLab.Solver_Manager.Model_Factory.Pre_Calculation.check_dependencies(Dict(1 => [1,
                                                                                            2]))
        return PeriLab.Data_Manager.get_properties(1, "Pre Calculation Model")
    end
    corr = Dict{String,Any}("Material Model" => "Correspondence Elastic",
                            "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                            "Shear Modulus" => 1.0)
    p = ut_dependencies(merge(corr, Dict{String,Any}("Bond Associated" => true)))
    @test p["Bond Associated Correspondence"] && p["Deformed Bond Geometry"]
    @test !haskey(p, "Shape Tensor")
    p = ut_dependencies(corr)
    @test p["Shape Tensor"] && p["Deformation Gradient"]
    p = ut_dependencies(Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    @test p["Deformed Bond Geometry"] && !haskey(p, "Shape Tensor")
end

@testset "strain Hooke matrix" begin
    ut_reset(2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Symmetry" => "isotropic plane strain",
                                  "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0))
    legacy = ut_legacy_dict(Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                             "Symmetry" => "isotropic plane strain",
                                             "Bulk Modulus" => 10.0,
                                             "Shear Modulus" => 4.0); dof = 2)
    @test Matrix(BBASIS.hooke_matrix(m, 2)) ≈
          Matrix(BBASIS.get_Hooke_matrix(legacy, legacy["Symmetry"], 2))
end


@testset "matrix-based correspondence init takes the block material" begin
    MB = PeriLab.Solver_Manager.Correspondence_matrix_based
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                                  "Shear Modulus" => 1.0, "Zero Energy Control" => "Global"))
    @test hasmethod(MB.init_model, Tuple{Vector{Int64},typeof(m),Int64})
    @test length(methods(MB.init_model)) == 1   # the Dict method is gone
end

@testset "dispatch without the material dict" begin
    MF = PeriLab.Solver_Manager.Model_Factory
    ut_reset(3; nnodes = 2)
    PeriLab.Data_Manager.set_block_id_list([1, 2])
    PeriLab.Data_Manager.init_properties()
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    PeriLab.Data_Manager.set_block_material(1, m)
    @test MF.has_block_model(1, "Material Model")
    @test !MF.has_block_model(2, "Material Model")
    @test MF.block_model_parameters(1, "Material Model") === m
    @test !MF.has_block_model(1, "Damage Model")
    @test_logs (:error, "Block 2 has no material model defined.") @test_throws PeriLab.PeriLabError begin
        BMAT.init_model([1, 2], 2)
    end
    @test hasmethod(BMAT.compute_model, Tuple{Vector{Int64},typeof(m),Int64,Float64,Float64})
end

@testset "2D symmetry check" begin
    ut_reset(2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                                  "Shear Modulus" => 1.0); dof = 2)
    @test_logs (:error,
                "Model definition is missing; plane stress or plane strain has to be defined for 2D") @test_throws PeriLab.PeriLabError begin
        BMAT.check_material_symmetry(m, 2)
    end
    ok = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                   "Symmetry" => "isotropic plane strain",
                                   "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0); dof = 2)
    @test isnothing(BMAT.check_material_symmetry(ok, 2))
end

