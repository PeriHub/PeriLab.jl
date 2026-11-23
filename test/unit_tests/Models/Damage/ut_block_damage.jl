# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const BDAM = PeriLab.Solver_Manager.Model_Factory.Damage

function ut_typed_damage(raw)
    ctx = PeriLab.ParameterSpec.ParseContext()
    wb = PeriLab.ParameterSpec.parse_model(:damage,
                                           Dict{String,Any}(string(k) => v
                                                            for (k, v) in raw),
                                           "test", ctx; name_key = "Damage Model")
    PeriLab.ParameterSpec.report!(ctx)
    return BDAM.block_damage(wb)
end

function ut_bind_allocations(tables)
    PeriLab.Data_Manager.bind_dependent_tables!(tables)
    return @allocated PeriLab.Data_Manager.bind_dependent_tables!(tables)
end

@testset "block damage binds its tables" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    file = joinpath(mktempdir(), "crit.txt")
    write(file, "header: Temperature Critical_Value\n0.0 1.0\n100.0 2.0\n")
    d = ut_typed_damage(Dict("Damage Model" => "Critical Energy", "Critical Value" => file))
    @test d isa BDAM.BlockDamage
    @test length(d.tables) == 1
    @test_logs (:error,
                "Field \"Temperature\" required by $(d.tables[1].source) does not exist or is not a per-node Vector{Float64}.") @test_throws PeriLab.PeriLabError PeriLab.Data_Manager.bind_dependent_tables!(d.tables)
    _, temperature = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    temperature .= [0.0, 100.0]
    PeriLab.Data_Manager.bind_dependent_tables!(d.tables)
    @test PeriLab.ParameterSpec.value(d.base.critical_value, 1) ≈ 1.0
    @test PeriLab.ParameterSpec.value(d.base.critical_value, 2) ≈ 2.0
    # called before every compute: binding must not allocate
    @test ut_bind_allocations(d.tables) == 0
    @test isempty(ut_typed_damage(Dict("Damage Model" => "Critical Stretch",
                                       "Critical Value" => 0.1)).tables)
end

@testset "block damage slot" begin
    PeriLab.Data_Manager.initialize_data()
    @test PeriLab.Data_Manager.get_block_models(1).damage === nothing
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1))
    PeriLab.Data_Manager.set_block_models(1, PeriLab.Data_Manager.BlockModels(damage = d))
    @test PeriLab.Data_Manager.get_block_models(1).damage === d
    @test "Block Models" in PeriLab.Data_Manager.CHECKPOINT_EXCLUDED_KEYS
end

@testset "typed interface critical values" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_block_id_list([2, 3, 1])
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 1.0,
                             "Interblock Damage" =>
                                 Dict("Interblock Critical Value 1_2" => 0.2,
                                      "Interblock Critical Value 2_3" => 0.3,
                                      "Interblock Critical Value 2_1" => 0.4,
                                      "Interblock Critical Value 4_1" => 0.9)))
    BDAM.init_interface_crit_values(d, 1)
    m = PeriLab.Data_Manager.get_crit_values_matrix()
    @test size(m) == (3, 3, 3)
    @test m[1, 2, 1] == 0.2 && m[2, 3, 1] == 0.3 && m[2, 1, 1] == 0.4
    @test m[1, 1, 1] == 1.0 && m[3, 3, 2] == 1.0          # filled with Critical Value
    plain = ut_typed_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 1.0))
    PeriLab.Data_Manager.initialize_data()
    BDAM.init_interface_crit_values(plain, 1)
    @test PeriLab.Data_Manager.get_crit_values_matrix() == fill(-1, (1, 1, 1))
end

@testset "typed anisotropic critical values" begin
    PeriLab.Data_Manager.initialize_data()
    d = ut_typed_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                             "Critical Value" => 1.0,
                             "Anisotropic Damage" => Dict("Critical Value X" => 1.0,
                                  "Critical Value Y" => 2.0)))
    BDAM.init_aniso_crit_values(d.base.anisotropic_damage, 1, 2)
    BDAM.init_aniso_crit_values(d.base.anisotropic_damage, 2, 3)
    z = ut_typed_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                             "Critical Value" => 1.0,
                             "Anisotropic Damage" => Dict("Critical Value X" => 1.0,
                                  "Critical Value Y" => 2.0,
                                  "Critical Value Z" => 3.0)))
    BDAM.init_aniso_crit_values(z.base.anisotropic_damage, 3, 3)
    aniso = PeriLab.Data_Manager.get_aniso_crit_values()
    @test aniso[1] == [1.0, 2.0]
    @test aniso[2] == [1.0, 2.0, 2.0]          # Z defaults to Y
    @test aniso[3] == [1.0, 2.0, 3.0]
end

function ut_stretch_setup(deformed)
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(2)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn .= 1
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist",
                                                                   Int64)
    nlist[1] = [2]
    nlist[2] = [1]
    PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64) .= 1
    PeriLab.Data_Manager.create_constant_node_scalar_field("Update", Bool)
    PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    len = PeriLab.Data_Manager.create_constant_bond_scalar_state("Bond Length", Float64)
    len[1] .= 1.0
    len[2] .= 1.0
    _, dlen = PeriLab.Data_Manager.create_bond_scalar_state("Deformed Bond Length", Float64)
    dlen[1] .= deformed[1]
    dlen[2] .= deformed[2]
    # bonds along x; the stretch only changes the length
    S = PeriLab.Data_Manager.SVector{2,Float64}
    bg = PeriLab.Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64, 2)
    _, dbg = PeriLab.Data_Manager.create_bond_vector_state("Deformed Bond Geometry",
                                                           Float64, 2)
    for (i, sign) in ((1, 1.0), (2, -1.0))
        bg[i][1] = S(sign, 0.0)
        dbg[i][1] = S(sign * deformed[i], 0.0)
    end
    return nothing
end

function ut_stretch_damage(raw; deformed = (1.2, 0.9))
    ut_stretch_setup(deformed)
    d = ut_typed_damage(raw)
    PeriLab.Data_Manager.set_block_models(1, PeriLab.Data_Manager.BlockModels(damage = d))
    BDAM.Critical_Stretch.init_model([1, 2], d.model, d, 1)
    BDAM.Critical_Stretch.compute_model([1, 2], d.model, d, 1, 0.0, 1.0)
    bd = PeriLab.Data_Manager.get_bond_damage("NP1")
    return (bd[1][1], bd[2][1])
end

@testset "critical stretch on the typed interface" begin
    # stretches: node 1 +0.2, node 2 -0.1
    @test ut_stretch_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 0.15)) == (0.0, 1.0)
    @test ut_stretch_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 0.05, "Compression" => "full")) ==
          (0.0, 0.0)
    @test ut_stretch_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 0.05)) == (0.0, 1.0)
end

@testset "anisotropic energy needs Anisotropic Damage" begin
    ut_stretch_setup((1.0, 1.0))
    PeriLab.Data_Manager.create_constant_node_scalar_field("Horizon", Float64) .= 1.0
    d = ut_typed_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                             "Critical Value" => 1.0))
    @test_logs (:error,
                "Critical Energy Anisotropic requires Anisotropic Damage.") @test_throws PeriLab.PeriLabError BDAM.Critical_Energy_Aniso.init_model([
                                                                                                                                                        1,
                                                                                                                                                        2
                                                                                                                                                    ],
                                                                                                                                                    d.model,
                                                                                                                                                    d,
                                                                                                                                                    1)
end

@testset "damage dispatch reads the block damage" begin
    ut_stretch_setup((1.2, 0.9))
    PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64) .= 1.0
    PeriLab.Data_Manager.create_node_scalar_field("Damage", Float64)
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch",
                             "Critical Value" => 0.15))
    PeriLab.Data_Manager.set_block_models(1, PeriLab.Data_Manager.BlockModels(damage = d))
    BDAM.init_model([1, 2], 1)
    BDAM.compute_model([1, 2], d, 1, 0.0, 1.0)
    @test PeriLab.Data_Manager.get_damage("NP1") == [1.0, 0.0]
    @test length(methods(BDAM.Critical_Stretch.compute_model)) == 1
end

# one bond 1 -> 2 (node 2 mirrored), undeformed (1, 0); bond energies of node 1 and 2
function ut_bond_energies(compression, deformed, force)
    S = PeriLab.Data_Manager.SVector{2,Float64}
    ut_stretch_setup((hypot(deformed...), hypot(deformed...)))
    dbg = PeriLab.Data_Manager.get_field("Deformed Bond Geometry", "NP1")
    dbg[1][1] = S(deformed)
    dbg[2][1] = -S(deformed)
    forces = PeriLab.Data_Manager.create_constant_bond_vector_state("Bond Forces", Float64,
                                                                    2)
    forces[1][1] = S(force)
    forces[2][1] = -S(force)
    PeriLab.Data_Manager.create_constant_node_scalar_field("Quad Horizon", Float64) .= 1.0
    inverse = PeriLab.Data_Manager.create_constant_bond_scalar_state("Inverse Bond ID",
                                                                     Int64)
    inverse[1] .= 1
    inverse[2] .= 1
    energy = PeriLab.Data_Manager.create_constant_bond_scalar_state("Bond Energy", Float64)
    PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64) .= 1.0
    PeriLab.Data_Manager.create_constant_node_scalar_field("Dissipated Energy", Float64)
    BDAM.Critical_Energy._critical_energy!([1, 2], compression == "none",
                                           compression == "shear", (i, j, k) -> Inf)
    return (energy[1][1], energy[2][1])
end

@testset "Compression of Critical Energy" begin
    modes = ("full", "none", "shear")
    # Δf = 2 force, Δu = deformed - (1, 0)
    expected(deformed, force) = 0.25 * abs(2 * (force[1] * (deformed[1] - 1) +
                                     force[2] * deformed[2]))
    # pure compression along the bond
    d, f = [0.9, 0.0], [-1.0, 0.0]
    e = Dict(m => ut_bond_energies(m, d, f) for m in modes)
    @test e["full"][1] ≈ expected(d, f) && e["full"][1] > 0
    @test e["full"][1] ≈ e["full"][2]
    @test e["none"] == (0.0, 0.0)
    @test all(abs.(e["shear"]) .< 1e-14)
    # pure sliding: Δu and Δf normal to the (slightly shortened) deformed bond
    θ = 0.1
    d, f = cos(θ) * [cos(θ), sin(θ)], [-sin(θ), cos(θ)]
    e = Dict(m => ut_bond_energies(m, d, f) for m in modes)
    @test e["full"][1] ≈ expected(d, f) && e["full"][1] > 0
    @test all(e["shear"] .≈ e["full"])
    @test e["none"] == (0.0, 0.0)
    # stretched bond: every mode counts it fully
    d, f = [1.1, 0.05], [1.0, 0.3]
    e = Dict(m => ut_bond_energies(m, d, f) for m in modes)
    @test e["full"][1] ≈ expected(d, f) && e["full"][1] > 0
    @test e["none"] == e["full"] && e["shear"] == e["full"]
end

@testset "Compression of Critical Stretch" begin
    damage(compression, deformed;
           crit = 0.05) = ut_stretch_damage(Dict("Damage Model" => "Critical Stretch",
                                                 "Critical Value" => crit,
                                                 "Compression" => compression);
                                            deformed = deformed)
    # node 1 stretched by 0.2, node 2 compressed along the bond by 0.1
    @test damage("full", (1.2, 0.9)) == (0.0, 0.0)
    @test damage("none", (1.2, 0.9)) == (0.0, 1.0)
    @test damage("shear", (1.2, 0.9)) == (0.0, 1.0)
    # sliding of a shortened bond: |Δu_t| / |ξ| = sin(0.3) ≈ 0.296 > 0.2 > stretch 0.045
    θ = 0.3
    slide(compression) = begin
        ut_stretch_setup((cos(θ), cos(θ)))
        dbg = PeriLab.Data_Manager.get_field("Deformed Bond Geometry", "NP1")
        dbg[1][1] = cos(θ) * PeriLab.Data_Manager.SVector(cos(θ), sin(θ))
        dbg[2][1] = -dbg[1][1]
        d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 0.2, "Compression" => compression))
        BDAM.Critical_Stretch.init_model([1, 2], d.model, d, 1)
        BDAM.Critical_Stretch.compute_model([1, 2], d.model, d, 1, 0.0, 1.0)
        bd = PeriLab.Data_Manager.get_bond_damage("NP1")
        (bd[1][1], bd[2][1])
    end
    @test slide("shear") == (0.0, 0.0)
    @test slide("none") == (1.0, 1.0)
    @test slide("full") == (1.0, 1.0)
end
