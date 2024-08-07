# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using TestSetExtensions
using Aqua
using Logging
using MPI
using PeriLab
Logging.disable_logging(Logging.Error)

# Aqua.test_all(PeriLab, ambiguities=false, stale_deps=(ignore=[:ZipArchives],))

include("helper.jl")

MPI.Init()

@testset ExtendedTestSet "PeriLab" begin

    @testset "unit_tests" begin
        @testset "ut_perilab" begin
            include("unit_tests/ut_perilab.jl")
        end
        @testset "FEM" begin
            @testset "ut_FEM_routines" begin
                include("unit_tests/FEM/ut_FEM_routines.jl")
            end

            @testset "ut_FEM_Factory" begin
                include("unit_tests/FEM/ut_FEM_Factory.jl")
            end
            @testset "ut_lagrange_element" begin
                include("unit_tests/FEM/Element_formulation/ut_lagrange_element.jl")
            end
        end
        @testset "Compute" begin
            @testset "ut_compute_global_values" begin
                include("unit_tests/Compute/ut_compute_global_values.jl")
            end
            @testset "ut_compute_field_values" begin
                include("unit_tests/Compute/ut_compute_field_values.jl")
            end
        end
        @testset "Support" begin
            @testset "Parameters" begin
                @testset "ut_parameter_handling" begin
                    include("unit_tests/Support/Parameters/ut_parameter_handling.jl")
                end
            end

            @testset "ut_helpers" begin
                include("unit_tests/Support/ut_helpers.jl")
            end

            @testset "ut_geometry" begin
                include("unit_tests/Support/ut_geometry.jl")
            end
        end
        @testset "Core" begin
            @testset "ut_data_manager" begin
                include("unit_tests/Core/ut_data_manager.jl")
            end
            @testset "ut_influence_function" begin
                include("unit_tests/Core/ut_Influence_function.jl")
            end
            @testset "Solver" begin

                @testset "ut_Solver_control" begin
                    include("unit_tests/Core/Solver/ut_Solver_control.jl")
                end

                @testset "ut_Verlet" begin
                    include("unit_tests/Core/Solver/ut_Verlet.jl")
                end

            end
            @testset "Module_inclusion" begin
                include("unit_tests/Core/Module_inclusion/ut_set_Modules.jl")
            end
            @testset "ut_BC_manager" begin
                include("unit_tests/Core/ut_BC_manager.jl")
            end
        end
        @testset "IO" begin

            @testset "ut_exodus_export" begin
                include("unit_tests/IO/ut_exodus_export.jl")
            end
            @testset "ut_read_inputdeck" begin
                include("unit_tests/IO/ut_read_inputdeck.jl")
            end
            @testset "ut_IO" begin
                include("unit_tests/IO/ut_IO.jl")
            end

            @testset "ut_mesh_data" begin
                include("unit_tests/IO/ut_mesh_data.jl")
            end

            @testset "ut_bond_filter" begin
                include("unit_tests/IO/ut_bond_filter.jl")
            end
            @testset "ut_logging" begin
                include("unit_tests/IO/ut_logging.jl")
            end
        end
        @testset "MPI" begin
            @testset "ut_MPI" begin
                include("unit_tests/MPI_communication/ut_MPI_call.jl")
            end

        end
        @testset "Physics" begin
            @testset "ut_templates" begin
                include("unit_tests/Physics/ut_templates.jl")
            end
            @testset "Additive" begin
                @testset "ut_Additive_Factory" begin
                    include("unit_tests/Physics/Additive/ut_Additive_Factory.jl")
                end
            end
            @testset "Corrosion" begin
                @testset "ut_Corrosion_Factory" begin
                    include("unit_tests/Physics/Corrosion/ut_Corrosion_Factory.jl")
                end
            end
            @testset "Thermal" begin
                @testset "ut_Thermal_Factory" begin
                    include("unit_tests/Physics/Thermal/ut_Thermal_Factory.jl")
                end
                @testset "ut_Thermal_Flow" begin
                    include("unit_tests/Physics/Thermal/ut_Thermal_flow.jl")
                end
                @testset "ut_Thermal_Expansion" begin
                    include("unit_tests/Physics/Thermal/ut_Thermal_expansion.jl")
                end
                @testset "ut_Heat_transfer" begin
                    include("unit_tests/Physics/Thermal/ut_Heat_transfer.jl")
                end
            end

            @testset "ut_Physics_Factory" begin
                include("unit_tests/Physics/ut_Physics_Factory.jl")
            end
            @testset "ut_Pre_calculation" begin
                include("unit_tests/Physics/Pre_calculation/ut_bond_deformation_gradient.jl")
            end
            @testset "ut_Damage" begin
                include("unit_tests/Physics/Damage/ut_Damage_Factory.jl")
                include("unit_tests/Physics/Damage/ut_Energy_release.jl")
            end
            @testset "ut_Material" begin
                @testset "ut_Material_Factory" begin
                    include("unit_tests/Physics/Material/ut_Material_Factory.jl")
                end
                @testset "ut_control" begin
                    include("unit_tests/Physics/Material/Zero_Energy_Control/ut_global_control.jl")
                end
                @testset "ut_material_basis" begin
                    include("unit_tests/Physics/Material/ut_material_basis.jl")
                end
                @testset "ut_bond_based" begin
                    include("unit_tests/Physics/Material/Material_Models/BondBased/ut_Bondbased_Elastic.jl")
                end
                @testset "ut_correspondence" begin
                    include("unit_tests/Physics/Material/Material_Models/ut_Correspondence.jl")
                    include("unit_tests/Physics/Material/Material_Models/ut_Correspondence_Plastic.jl")
                    include("unit_tests/Physics/Material/Material_Models/ut_Correspondence_UMAT.jl")
                    include("unit_tests/Physics/Material/Material_Models/ut_Bond_Associated_Correspondence.jl")
                end
                @testset "ut_ordinary" begin
                    include("unit_tests/Physics/Material/Material_Models/Ordinary/ut_ordinary.jl")
                    include("unit_tests/Physics/Material/Material_Models/ut_PD_Solid_Elastic.jl")
                    include("unit_tests/Physics/Material/Material_Models/ut_PD_Solid_Plastic.jl")
                end
            end
        end
    end

    @testset "fullscale_tests" begin
        @testset "test_reload" begin
            include("fullscale_tests/test_reload/test_reload.jl")
        end
        @testset "test_dry_run" begin
            include("fullscale_tests/test_dry_run/test_dry_run.jl")
        end
        @testset "test_additive_simple" begin
            include("fullscale_tests/test_additive/test_additive.jl")
        end
        @testset "test_test_bond_based_elastic" begin
            include("fullscale_tests/test_bond_based_elastic/test_bond_based_elastic.jl")
        end
        @testset "test_heat_transfer" begin
            include("fullscale_tests/test_heat_transfer/test_heat_transfer.jl")
        end
        @testset "test_BCs" begin
            include("fullscale_tests/test_BCs/test_BCs.jl")
        end
        # @testset "test_body_force" begin
        #     include("fullscale_tests/test_body_force/test_body_force.jl")
        # end
        @testset "test_contact" begin
            include("fullscale_tests/test_contact/test_contact.jl")
        end
        @testset "test_PD_Solid_Elastic" begin
            include("fullscale_tests/test_PD_solid_elastic/test_PD_solid_elastic.jl")
        end
        @testset "test_PD_Solid_Elastic_3D" begin
            include("fullscale_tests/test_PD_solid_elastic_3D/test_PD_solid_elastic_3D.jl")
        end
        # @testset "test_PD_Solid_Plastic" begin
        #     include("fullscale_tests/test_PD_solid_plastic/test_PD_solid_plastic.jl")
        # end
        @testset "test_Critical_stretch" begin
            include("fullscale_tests/test_critical_stretch/test_critical_stretch.jl")
        end
        @testset "test_critical_energy" begin
            include("fullscale_tests/test_critical_energy/test_critical_energy.jl")
        end
        @testset "test_thermal_expansion" begin
            include("fullscale_tests/test_thermal_expansion/test_thermal_expansion.jl")
        end
        @testset "test_thermal_flow" begin
            include("fullscale_tests/test_thermal_flow/test_thermal_flow.jl")
            # include("fullscale_tests/test_thermal_flow_paper/test_thermal_flow.jl")
        end
        @testset "test_Correspondence_Elastic" begin
            include("fullscale_tests/test_correspondence_elastic/test_correspondence_elastic.jl")
        end
        @testset "test_Correspondence_Elastic_Plastic" begin
            include("fullscale_tests/test_correspondence_elastic_plastic/test_correspondence_elastic_plastic.jl")
        end
        @testset "test_Correspondence_Elastic_with_zero_E_control" begin
            include("fullscale_tests/test_correspondence_elastic_with_zero_E_control/test_correspondence_elastic_with_zero_E_control.jl")
        end
        @testset "test_Correspondence_Elastic" begin
            include("fullscale_tests/test_correspondence_elastic_3D/test_correspondence_elastic_3D.jl")
        end
        @testset "test_DCB" begin
            include("fullscale_tests/test_DCB/test_DCB.jl")
        end
        @testset "test_Abaqus" begin
            include("fullscale_tests/test_Abaqus/test_Abaqus.jl")
        end
        @testset "test_aniso_damage" begin
            include("fullscale_tests/test_aniso_damage/test_aniso_damage.jl")
        end
        @testset "test_material_field" begin
            include("fullscale_tests/test_material_field/test_material_field.jl")
        end
        @testset "test_FEM" begin
            include("fullscale_tests/test_FEM/test_FEM.jl")
        end
    end


end

MPI.Finalize()

#cleanup
rm("tmp", force=true, recursive=true)
run(`find . -name "*.log" -type f -delete`)