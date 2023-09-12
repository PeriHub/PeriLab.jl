# include("../src/PeriLab.jl")
# import .PeriLab
using Test
using Aqua
using Logging
using PeriLab
Logging.disable_logging(Logging.Error)
MPI.Init()

@testset "Support" begin

    @testset "Parameters" begin

        @testset "ut_parameter_handling" begin
            include("../src/Support/Parameters/unit_test/ut_parameter_handling.jl")
        end

    end

    @testset "ut_data_manager" begin
        include("../src/Support/unit_test/ut_data_manager.jl")
    end

    @testset "ut_helpers" begin
        include("../src/Support/unit_test/ut_helpers.jl")
    end

    @testset "ut_tools" begin
        include("../src/Support/unit_test/ut_tools.jl")
    end

    @testset "ut_geometry" begin
        include("../src/Support/unit_test/ut_geometry.jl")
    end
end

@testset "Core" begin


    @testset "Solver" begin

        # @testset "ut_Solver_control" begin
        #     include("../src/Core/Solver/unit_test/ut_Solver_control.jl")
        # end

        # @testset "ut_Verlet" begin
        # include("../src/Core/Solver/unit_test/ut_Verlet.jl")
        # end

    end

    @testset "ut_BC_manager" begin
        include("../src/Core/unit_test/ut_BC_manager.jl")
    end
end

@testset "IO" begin

    @testset "ut_exodus_export" begin
        include("../src/IO/unit_test/ut_exodus_export.jl")
    end

    @testset "ut_IO" begin
        include("../src/IO/unit_test/ut_IO.jl")
    end

    @testset "ut_mesh_data" begin
        include("../src/IO/unit_test/ut_mesh_data.jl")
    end
end


@testset "MPI" begin

    @testset "ut_MPI" begin
        #   include("../src/MPI_communication/unit_test/ut_MPI.jl")
    end

end


@testset "Physics" begin

    @testset "Thermal" begin end
    @testset "Correspondence" begin


        include("../src/Physics/Thermal/Correspondence/unit_test/ut_Thermal_correspondence.jl")



    end

    @testset "ut_Physics_Factory" begin
        include("../src/Physics/unit_test/ut_Physics_Factory.jl")
    end

end

MPI.Finalize()
# Aqua.test_all(PeriLab)
