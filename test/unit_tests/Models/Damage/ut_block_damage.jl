# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const BDAM = PeriLab.Solver_Manager.Model_Factory.Damage

function ut_typed_damage(raw)
    ctx = PeriLab.ParameterSpec.ParseContext()
    wb = PeriLab.ParameterSpec.parse_model(:damage,
                                           Dict{String,Any}(string(k) => v for (k, v) in raw),
                                           "test", ctx; name_key = "Damage Model")
    PeriLab.ParameterSpec.report!(ctx)
    return BDAM.block_damage(wb)
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
    @test isempty(ut_typed_damage(Dict("Damage Model" => "Critical Stretch",
                                       "Critical Value" => 0.1)).tables)
end

@testset "block damage slot" begin
    PeriLab.Data_Manager.initialize_data()
    @test PeriLab.Data_Manager.get_block_damage(1) === nothing
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1))
    PeriLab.Data_Manager.set_block_damage(1, d)
    @test PeriLab.Data_Manager.get_block_damage(1) === d
    @test "Block Damages" in PeriLab.Data_Manager.CHECKPOINT_EXCLUDED_KEYS
end

@testset "typed interface critical values" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_block_id_list([2, 3, 1])
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 1.0,
                             "Interblock Damage" => Dict("Interblock Critical Value 1_2" => 0.2,
                                                         "Interblock Critical Value 2_3" => 0.3,
                                                         "Interblock Critical Value 2_1" => 0.4,
                                                         "Interblock Critical Value 4_1" => 0.9)))
    BDAM.init_interface_crit_values(d, 1)
    m = PeriLab.Data_Manager.get_crit_values_matrix()
    @test size(m) == (3, 3, 3)
    @test m[1, 2, 1] == 0.2 && m[2, 3, 1] == 0.3 && m[2, 1, 1] == 0.4
    @test m[1, 1, 1] == 1.0 && m[3, 3, 2] == 1.0          # filled with Critical Value
    plain = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 1.0))
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
