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
MPI.Init()

@testset ExtendedTestSet "PeriLab" begin

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

        @testset "ut_tools" begin
            @includetests["unit_tests/Support/ut_tools"]
        end

        @testset "ut_geometry" begin
            @includetests["unit_tests/Support/ut_geometry"]
        end
    end

    @testset "Core" begin


        @testset "Solver" begin

            # @testset "ut_Solver_control" begin
            #     include("unit_tests/Core/Solver//ut_Solver_control.jl")
            # end

            # @testset "ut_Verlet" begin
            # include("unit_tests/Core/Solver//ut_Verlet.jl")
            # end

        end

        @testset "ut_BC_manager" begin
            @includetests["unit_tests/Core/ut_BC_manager"]
        end
    end

    @testset "IO" begin

        @testset "ut_exodus_export" begin
            @includetests["unit_tests/IO/ut_exodus_export"]
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
    end
    @testset "Compute" begin
        @includetests["unit_tests/Compute/ut_compute_force"]
    end

    @testset "MPI" begin

        @testset "ut_MPI" begin
            #   include("unit_tests/MPI_communication/ut_MPI.jl")
        end

    end


    @testset "Physics" begin

        @testset "Thermal" begin end

        @testset "ut_Physics_Factory" begin
            @includetests["unit_tests/Physics/ut_Physics_Factory"]
        end
        @testset "ut_Damage" begin
            @includetests["unit_tests/Physics/Damage/ut_Damage_Factory"]
        end
        @testset "ut_Material" begin
            @testset "ut_material_basis" begin
                @includetests["unit_tests/Physics/Material/ut_material_basis"]
            end
            @testset "ut_correspondence" begin
                @includetests["unit_tests/Physics/Material/Material_Models/ut_Correspondence"]
            end
            @testset "ut_ordinary" begin
                @includetests["unit_tests/Physics/Material/Material_Models/Ordinary/ut_ordinary"]
            end
        end
    end

    @testset "fullscale_tests" begin
        @testset "test_BCs" begin
            @includetests["fullscale_tests/test_BCs/test_BCs"]
        end
        @testset "test_PD_Solid_Elastic" begin
            @includetests["fullscale_tests/test_PD_Solid_Elastic/test_PD_Solid_Elastic"]
        end

        @testset "test_Critical_stretch" begin
            @includetests["fullscale_tests/test_Critical_stretch/test_Critical_stretch"]
        end
    end

    # @testset "test_Correspondence_Elastic" begin
    #     @includetests["fullscale_tests/test_Correspondence_Elastic/test_Correspondence_Elastic"]
    # end
end

MPI.Finalize()
# Aqua.test_all(PeriLab)
