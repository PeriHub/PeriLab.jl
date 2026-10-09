# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# The templates are what developers copy to write a new model: their functions
# must take exactly the arguments the factory passes to a real model of the
# same category.

const UT_MF = PeriLab.Solver_Manager.Model_Factory

# argument types of the single method of `f`, with the model's own parameter
# struct `P` replaced by :params
function ut_signature(f, P)
    ms = collect(methods(f))
    @test length(ms) == 1
    return [T === P ? :params : T for T in ms[1].sig.parameters[2:end]]
end

# the template accepts every call the real model gets, function by function
function ut_same_interface(template, P_template, model, P_model, functions)
    for f in functions
        @test isdefined(template, f)
        t = ut_signature(getfield(template, f), P_template)
        m = ut_signature(getfield(model, f), P_model)
        @test length(t) == length(m) &&
              all(a === :params ? b === :params : b !== :params && b <: a for (a, b) in zip(t, m))
    end
end

const UT_NODES = Vector{Int64}(1:3)

PeriLab.Data_Manager.initialize_data()
PeriLab.Data_Manager.set_num_controller(3)
PeriLab.Data_Manager.set_dof(2)

@testset "ut_material_template" begin
    MT = UT_MF.Material.Material_template
    BE = UT_MF.Material.Bondbased_Elastic
    ut_same_interface(MT, MT.MaterialTemplateParams, BE, BE.BondbasedElasticParams,
                      (:init_model, :compute_model, :fields_for_local_synchronization))
    @test !MT.fe_support()
    p = MT.MaterialTemplateParams()
    material = typed_block_material(Dict("Material Model" => "Bond-based Elastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    MT.init_model(UT_NODES, p, material)
    MT.compute_model(UT_NODES, p, material, 1, 0.0, 0.0)
    MT.fields_for_local_synchronization("Material Model")
end

@testset "ut_correspondence_template" begin
    # not loaded by PeriLab (it lives next to the material template); load it
    # where a copy would live, beside the correspondence models
    C = UT_MF.Material.Correspondence
    isdefined(C, :Correspondence_template) ||
        Base.include(C,
                     joinpath(pkgdir(PeriLab), "src", "Models", "Material", "Material_Models",
                              "Material_template", "correspondence_template.jl"))
    CT = C.Correspondence_template
    CE = C.Correspondence_Elastic
    ut_same_interface(CT, CT.CorrespondenceTemplateParams, CE, CE.CorrespondenceElasticParams,
                      (:init_model, :compute_stresses, :compute_stresses_ba,
                       :fields_for_local_synchronization))
    @test !CT.fe_support()
    p = CT.CorrespondenceTemplateParams()
    material = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    CT.init_model(UT_NODES, p, material)
    stress_NP1 = zeros(1, 2, 2)
    stress_NP1[1, :, 1] = [-1, 2.2]
    vec = CT.compute_stresses(Vector{Int64}(1:1), 2, p, material, 0.0, 0.0, ones(1, 2, 2),
                              zeros(1, 2, 2), stress_NP1)
    @test vec === stress_NP1
    @test vec[1, :, 1] == [-1, 2.2]
    @test_logs (:error,
                "Correspondence Template not yet implemented for bond associated.") @test_throws PeriLab.PeriLabError CT.compute_stresses_ba(Vector{Int64}(1:1),
                                                                                                                                             [[1]],
                                                                                                                                             2,
                                                                                                                                             p,
                                                                                                                                             material,
                                                                                                                                             0.0,
                                                                                                                                             0.0,
                                                                                                                                             nothing,
                                                                                                                                             nothing,
                                                                                                                                             nothing)
    CT.fields_for_local_synchronization("Material Model")
end

@testset "ut_damage_template" begin
    DT = UT_MF.Damage.Damage_template
    CS = UT_MF.Damage.Critical_Stretch
    ut_same_interface(DT, DT.DamageTemplateParams, CS, CS.CriticalStretchParams,
                      (:init_model, :compute_model, :fields_for_local_synchronization))
    p = DT.DamageTemplateParams()
    damage = UT_MF.Damage.BlockDamage(nothing, p, PeriLab.ParameterSpec.Table1D[])
    DT.init_model(UT_NODES, p, damage, 1)
    DT.compute_model(UT_NODES, p, damage, 1, 0.0, 0.0)
    DT.fields_for_local_synchronization("Damage Model")
end

@testset "ut_thermal_template" begin
    TT = UT_MF.Thermal.Thermal_template
    HT = UT_MF.Thermal.Heat_Transfer
    ut_same_interface(TT, TT.ThermalTemplateParams, HT, HT.HeatTransferParams,
                      (:init_model, :compute_model, :fields_for_local_synchronization))
    p = TT.ThermalTemplateParams()
    TT.init_model(UT_NODES, p, nothing, 1)
    TT.compute_model(UT_NODES, p, nothing, 1, 0.0, 0.0)
    TT.fields_for_local_synchronization("Thermal Model")
end

@testset "ut_additive_template" begin
    # no additive model ships with PeriLab (they are licensed); the factory calls
    # them like degradation models
    AT = UT_MF.Additive.Additive_template
    BC = UT_MF.Degradation.Bondbased_Corrosion
    ut_same_interface(AT, AT.AdditiveTemplateParams, BC, BC.BondbasedCorrosionParams,
                      (:init_model, :compute_model, :fields_for_local_synchronization))
    p = AT.AdditiveTemplateParams()
    AT.init_model(UT_NODES, p, 1)
    AT.compute_model(UT_NODES, p, 1, 0.0, 0.0)
    AT.fields_for_local_synchronization("Additive Model")
end

@testset "ut_degradation_template" begin
    GT = UT_MF.Degradation.Degradation_template
    BC = UT_MF.Degradation.Bondbased_Corrosion
    ut_same_interface(GT, GT.DegradationTemplateParams, BC, BC.BondbasedCorrosionParams,
                      (:init_model, :compute_model, :fields_for_local_synchronization))
    p = GT.DegradationTemplateParams()
    GT.init_model(UT_NODES, p, 1)
    GT.compute_model(UT_NODES, p, 1, 0.0, 0.0)
    GT.fields_for_local_synchronization("Degradation Model")
end

@testset "ut_pre_calculation_template" begin
    PT = UT_MF.Pre_Calculation.Pre_calculation_template
    ST = UT_MF.Pre_Calculation.Shape_Tensor
    ut_same_interface(PT, nothing, ST, nothing,
                      (:init_model, :compute, :fields_for_local_synchronization))
    PT.init_model(UT_NODES, 1)
    PT.compute(UT_NODES, 1)
    PT.fields_for_local_synchronization("Pre Calculation Model")
end

@testset "ut_contact_template" begin
    CT = UT_MF.Contact.Contact_template
    PM = UT_MF.Contact.Penalty_Model
    ut_same_interface(CT, CT.ContactTemplateParams, PM, PM.PenaltyContactParams,
                      (:contact_model_name, :init_contact_model, :compute_contact_model))
    p = CT.ContactTemplateParams()
    CT.init_contact_model(p, nothing)
    CT.compute_contact_model("cg", p, nothing, (m, s, f) -> nothing, (s, m, f) -> nothing)
end

@testset "ut_FEM_template" begin
    FT = PeriLab.Solver_Manager.FEM.FEM_template
    LE = PeriLab.Solver_Manager.FEM.Lagrange_element
    ut_same_interface(FT, nothing, LE, nothing,
                      (:element_name, :init_element, :create_element_matrices))
    FT.init_element(UT_NODES, nothing, [1, 1])
    @test_throws PeriLab.PeriLabError FT.create_element_matrices(2, [2, 2], [1, 1],
                                                                 zeros(2, 2), zeros(2, 2))
end
