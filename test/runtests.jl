# include("../src/PeriLab.jl")
# import .PeriLab
using Test
using TestSetExtensions
using Aqua
using Logging
using PeriLab
using MPI
Logging.disable_logging(Logging.Error)
MPI.Init()

@testset ExtendedTestSet "PeriLab" begin

    @testset "Support" begin

        @testset "Parameters" begin

            @testset "ut_parameter_handling" begin
                include("unit_tests/Support/Parameters/ut_parameter_handling.jl")
            end

        end

        @testset "ut_data_manager" begin
            include("unit_tests/Support/ut_data_manager.jl")
        end

        @testset "ut_helpers" begin
            include("unit_tests/Support/ut_helpers.jl")
        end

        @testset "ut_tools" begin
            include("unit_tests/Support/ut_tools.jl")
        end

        @testset "ut_geometry" begin
            include("unit_tests/Support/ut_geometry.jl")
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
            include("unit_tests/Core/ut_BC_manager.jl")
        end
    end

    @testset "IO" begin

        @testset "ut_exodus_export" begin
            include("unit_tests/IO/ut_exodus_export.jl")
        end

        @testset "ut_IO" begin
            include("unit_tests/IO/ut_IO.jl")
        end

        @testset "ut_mesh_data" begin
            include("unit_tests/IO/ut_mesh_data.jl")
        end
    end


    @testset "MPI" begin

        @testset "ut_MPI" begin
            #   include("unit_tests/MPI_communication/ut_MPI.jl")
        end

    end


    @testset "Physics" begin

        @testset "Thermal" begin end
        @testset "Correspondence" begin


            include("unit_tests/Physics/Thermal/Correspondence/ut_Thermal_correspondence.jl")



        end

        @testset "ut_Physics_Factory" begin
            include("unit_tests/Physics/ut_Physics_Factory.jl")
        end
        @testset "ut_Damage" begin
            include("unit_tests/Physics/Damage/ut_Damage_Factory.jl")
        end
    end

    @testset "Test_PD_Solid_Elastic" begin

        include("test_PD_Solid_Elastic/test_PD_Solid_Elastic.jl")

    end
end

MPI.Finalize()
# Aqua.test_all(PeriLab)
