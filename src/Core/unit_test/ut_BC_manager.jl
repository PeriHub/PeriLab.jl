include("../BC_manager.jl")
using Test

using .Boundary_conditions
@testset "ut_evaluation" begin
    unit = ones(Float32, 3)
    dof = 2
    time = 2
    coor = zeros(3, 3)
    bc = "10"
    @test (10 * unit) == Boundary_conditions.eval_bc(bc, coor, time, dof)
    bc = "x"
    for i in 1:4
        coor[3, 1] = i * i - 2
        @test (coor[:, 1]) == Boundary_conditions.eval_bc(bc, coor, time, dof)
    end
    dof = 3
    bc = "t"
    @test (time * unit) == Boundary_conditions.eval_bc(bc, coor, time, dof)
    bc = "t*x"
    @test (time .* coor[:, 1]) == Boundary_conditions.eval_bc(bc, coor, time, dof)
    for t in 0:4
        bc = "if t>2 0 else 20 end"
        if t > 2
            @test (0.0 * unit) == Boundary_conditions.eval_bc(bc, coor, t, dof)
        else
            @test (20.0 * unit) == Boundary_conditions.eval_bc(bc, coor, t, dof)
        end
    end
end