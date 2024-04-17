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

Aqua.test_all(PeriLab, ambiguities=false, stale_deps=(ignore=[:ZipArchives],))

include("helper.jl")

MPI.Init()

@testset ExtendedTestSet "PeriLab" begin

    @testset "unit_tests" begin
        @testset "ut_perilab" begin
            @includetests["unit_tests/ut_perilab"]
        end
        @testset "FEM" begin
            @testset "ut_FEM_routines" begin
                @includetests["unit_tests/FEM/ut_FEM_routines"]
            end

            @testset "ut_FEM_Factory" begin
                @includetests["unit_tests/FEM/ut_FEM_Factory"]
            end
            @testset "ut_lagrange_element" begin
                @includetests["unit_tests/FEM/Element_formulation/ut_lagrange_element"]
            end
        end
        @testset "Compute" begin
            @testset "ut_compute_global_values" begin
                @includetests["unit_tests/Compute/ut_compute_global_values"]
            end
            @testset "ut_compute_field_values" begin
                @includetests["unit_tests/Compute/ut_compute_field_values"]
            end
        end
        @testset "Support" begin
            @testset "Parameters" begin
                @testset "ut_parameter_handling" begin
                    @includetests["unit_tests/Support/Parameters/ut_parameter_handling"]
                end
            end

            @testset "ut_data_manager" begin
                @includetests["unit_tests/Support/ut_data_manager"]
            end

            @testset "ut_helpers" begin
                @includetests["unit_tests/Support/ut_helpers"]
            end

            @testset "ut_geometry" begin
                @includetests["unit_tests/Support/ut_geometry"]
            end
        end
        @testset "Core" begin
            @testset "Solver" begin

                @testset "ut_Solver_control" begin
                    @includetests["unit_tests/Core/Solver/ut_Solver_control"]
                end

                @testset "ut_Verlet" begin
                    @includetests["unit_tests/Core/Solver/ut_Verlet"]
                end

            end
            @testset "Module_inclusion" begin
                @includetests["unit_tests/Core/Module_inclusion/ut_set_Modules"]
            end
            @testset "ut_BC_manager" begin
                @includetests["unit_tests/Core/ut_BC_manager"]
            end
        end
        @testset "IO" begin

            @testset "ut_exodus_export" begin
                @includetests["unit_tests/IO/ut_exodus_export"]
            end
            @testset "ut_read_inputdeck" begin
                @includetests["unit_tests/IO/ut_read_inputdeck"]
            end
            @testset "ut_IO" begin
                @includetests["unit_tests/IO/ut_IO"]
            end

            @testset "ut_mesh_data" begin
                @includetests["unit_tests/IO/ut_mesh_data"]
            end

            @testset "ut_bond_filter" begin
                @includetests["unit_tests/IO/ut_bond_filter"]
            end
            @testset "ut_logging" begin
                @includetests["unit_tests/IO/ut_logging"]
            end
        end
        @testset "MPI" begin
            @testset "ut_MPI" begin
                @includetests["unit_tests/MPI_communication/ut_MPI_call"]
            end

        end
        @testset "Physics" begin
            @testset "ut_templates" begin
                @includetests["unit_tests/Physics/ut_templates"]
            end
            @testset "Additive" begin
                @testset "ut_Additive_Factory" begin
                    @includetests["unit_tests/Physics/Additive/ut_Additive_Factory"]
                end
            end
            @testset "Corrosion" begin
                @testset "ut_Corrosion_Factory" begin
                    @includetests["unit_tests/Physics/Corrosion/ut_Corrosion_Factory"]
                end
            end
            @testset "Thermal" begin
                @testset "ut_Thermal_Factory" begin
                    @includetests["unit_tests/Physics/Thermal/ut_Thermal_Factory"]
                end
                @testset "ut_Thermal_Flow" begin
                    @includetests["unit_tests/Physics/Thermal/ut_Thermal_flow"]
                end
                @testset "ut_Thermal_Expansion" begin
                    @includetests["unit_tests/Physics/Thermal/ut_Thermal_expansion"]
                end
                @testset "ut_Heat_transfer" begin
                    @includetests["unit_tests/Physics/Thermal/ut_Heat_transfer"]
                end
            end

            @testset "ut_Physics_Factory" begin
                @includetests["unit_tests/Physics/ut_Physics_Factory"]
            end
            @testset "ut_Damage" begin
                @includetests["unit_tests/Physics/Damage/ut_Damage_Factory"]
                @includetests["unit_tests/Physics/Damage/ut_Energy_release"]
            end
            @testset "ut_Material" begin
                @testset "ut_Material_Factory" begin
                    @includetests["unit_tests/Physics/Material/ut_Material_Factory"]
                end
                @testset "ut_control" begin
                    @includetests["unit_tests/Physics/Material/Zero_Energy_Control/ut_global_control"]
                end
                @testset "ut_material_basis" begin
                    @includetests["unit_tests/Physics/Material/ut_material_basis"]
                end
                @testset "ut_correspondence" begin
                    @includetests["unit_tests/Physics/Material/Material_Models/ut_Correspondence"]
                    @includetests["unit_tests/Physics/Material/Material_Models/ut_Correspondence_UMAT"]
                end
                @testset "ut_ordinary" begin
                    @includetests["unit_tests/Physics/Material/Material_Models/Ordinary/ut_ordinary"]
                    @includetests["unit_tests/Physics/Material/Material_Models/ut_PD_Solid_Plastic"]
                end
            end
        end
    end

    @testset "fullscale_tests" begin
        @testset "test_reload" begin
            @includetests["fullscale_tests/test_reload/test_reload"]
        end
        @testset "test_dry_run" begin
            @includetests["fullscale_tests/test_dry_run/test_dry_run"]
        end
        @testset "test_additive_simple" begin
            @includetests["fullscale_tests/test_additive/test_additive"]
        end
        @testset "test_test_bond_based_elastic" begin
            @includetests["fullscale_tests/test_bond_based_elastic/test_bond_based_elastic"]
        end
        @testset "test_heat_transfer" begin
            @includetests["fullscale_tests/test_heat_transfer/test_heat_transfer"]
        end
        @testset "test_BCs" begin
            @includetests["fullscale_tests/test_BCs/test_BCs"]
        end
        # @testset "test_body_force" begin
        #     @includetests["fullscale_tests/test_body_force/test_body_force"]
        # end
        @testset "test_contact" begin
            @includetests["fullscale_tests/test_contact/test_contact"]
        end
        @testset "test_PD_Solid_Elastic" begin
            @includetests["fullscale_tests/test_PD_solid_elastic/test_PD_solid_elastic"]
        end
        @testset "test_PD_Solid_Elastic_3D" begin
            @includetests["fullscale_tests/test_PD_solid_elastic_3D/test_PD_solid_elastic_3D"]
        end
        @testset "test_Critical_stretch" begin
            @includetests["fullscale_tests/test_critical_stretch/test_critical_stretch"]
        end
        @testset "test_critical_energy" begin
            @includetests["fullscale_tests/test_critical_energy/test_critical_energy"]
        end
        @testset "test_thermal_expansion" begin
            @includetests["fullscale_tests/test_thermal_expansion/test_thermal_expansion"]
        end
        @testset "test_thermal_flow" begin
            @includetests["fullscale_tests/test_thermal_flow/test_thermal_flow"]
            # @includetests["fullscale_tests/test_thermal_flow_paper/test_thermal_flow"]
        end
        @testset "test_Correspondence_Elastic" begin
            @includetests["fullscale_tests/test_correspondence_elastic/test_correspondence_elastic"]
        end
        @testset "test_Correspondence_Elastic_Plastic" begin
            @includetests["fullscale_tests/test_correspondence_elastic_plastic/test_correspondence_elastic_plastic"]
        end
        @testset "test_Correspondence_Elastic_with_zero_E_control" begin
            @includetests["fullscale_tests/test_correspondence_elastic_with_zero_E_control/test_correspondence_elastic_with_zero_E_control"]
        end
        @testset "test_Correspondence_Elastic" begin
            @includetests["fullscale_tests/test_correspondence_elastic_3D/test_correspondence_elastic_3D"]
        end
        @testset "test_DCB" begin
            @includetests["fullscale_tests/test_DCB/test_DCB"]
        end
        @testset "test_Abaqus" begin
            @includetests["fullscale_tests/test_Abaqus/test_Abaqus"]
        end
        @testset "test_aniso_damage" begin
            @includetests["fullscale_tests/test_aniso_damage/test_aniso_damage"]
        end
        @testset "test_material_field" begin
            @includetests["fullscale_tests/test_material_field/test_material_field"]
        end
        @testset "test_FEM" begin
            @includetests["fullscale_tests/test_FEM/test_FEM"]
        end
    end


end

MPI.Finalize()

#cleanup
rm("tmp", force=true, recursive=true)
files = readdir()
log_files = filter(endswith(".log"), files)
for file in log_files
    rm(file)
end