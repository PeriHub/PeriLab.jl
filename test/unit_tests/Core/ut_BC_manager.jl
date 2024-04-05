# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../../../src/Core/BC_manager.jl")
# include("../../../src/Support/data_manager.jl")
# include("../../../src/Support/Parameters/parameter_handling.jl")
# using Reexport
# @reexport using .Parameter_Handling
# @reexport using .Boundary_conditions

using Test


@testset "ut_clean_up" begin
    @test PeriLab.Solver.Boundary_conditions.clean_up("") == ""
    @test PeriLab.Solver.Boundary_conditions.clean_up("-") == " .- "
    @test PeriLab.Solver.Boundary_conditions.clean_up("+") == " .+ "
    @test PeriLab.Solver.Boundary_conditions.clean_up("*") == " .* "
    @test PeriLab.Solver.Boundary_conditions.clean_up("/") == " ./ "
    @test PeriLab.Solver.Boundary_conditions.clean_up("5e-8") == "5e-8"
    @test PeriLab.Solver.Boundary_conditions.clean_up("5.0e-8") == "5.0e-8"
    @test PeriLab.Solver.Boundary_conditions.clean_up("5.0e-8-5") == "5.0e-8 .- 5"
    @test PeriLab.Solver.Boundary_conditions.clean_up("5.0e-8/5e+5*t") == "5.0e-8 ./ 5e+5 .* t"
    @test PeriLab.Solver.Boundary_conditions.clean_up("5-8.0/5e+5") == "5 .- 8.0 ./ 5e+5"
    @test PeriLab.Solver.Boundary_conditions.clean_up("2") == "2"
    @test PeriLab.Solver.Boundary_conditions.clean_up("2.") == "2."
    @test PeriLab.Solver.Boundary_conditions.clean_up("2.0") == "2.0"
    @test PeriLab.Solver.Boundary_conditions.clean_up("2.0.-1") == "2.0 .- 1"
end

@testset "ut_evaluation" begin
    unit = ones(Float64, 3)
    dof = 2
    time::Float64 = 2
    coor = zeros(3, 3)
    bc = Int64(10)
    values = ones(Float64, 3)
    @test (10 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, time, dof, false)
    bc = Float64(10)
    @test (10 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, time, dof, false)
    bc = Float64(10)
    @test (10 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, time, dof, false)
    bc = string(10)
    @test (10 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, time, dof, false)
    bc = "x"
    for i in 1:4
        coor[3, 1] = i * i - 2
        @test (coor[:, 1]) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, time, dof, false)
    end
    dof = 3
    bc = "t"
    @test (time * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, time, dof, false)
    bc = "t*x"
    @test (time .* coor[:, 1]) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, time, dof, false)
    for t in 0:4
        bc = "if t>2 0 else 20 end"
        if t > 2
            @test (0.0 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, Float64(t), dof, false)
        else
            @test (20.0 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, Float64(t), dof, false)
        end
    end
    for t in 0:2
        bc = "100.0"
        if t == 0
            @test (100.0 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, Float64(t), dof, false)
            @test (100.0 * unit) == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, Float64(t), dof, true)
        else
            @test unit == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, coor, Float64(t), dof, true)
        end
    end
    @test values == PeriLab.Solver.Boundary_conditions.eval_bc(values, bc, Matrix{Float64}(undef, 0, 0), time, dof, false)
end

@testset "ut_boundary_condition" begin
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.clear_data_manager()
    test_Data_manager.set_dof() = 2
    params = Dict()
    bcs = PeriLab.Solver.Boundary_conditions.boundary_condition(params, test_Data_manager)
    @test length(bcs) == 0
    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Variable" => "Forces", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Variable" => "Displacements", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "0"), "BC_3" => Dict("Variable" => "Displacements", "Node Set" => "Nset_3", "Coordinate" => "z", "Value" => "0")))

    test_Data_manager.set_nset("Nset_1", [1, 2, 3])
    test_Data_manager.set_nset("Nset_2", [3, 4, 7, 10])
    test_Data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 3, 3 => 4, 4 => 2, 5 => 5, 6 => 6, 7 => 7, 8 => 8, 9 => 9, 10 => 10))

    bcs = PeriLab.Solver.Boundary_conditions.boundary_condition(params, test_Data_manager)
    @test isnothing(bcs)
    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Variable" => "Forces", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Variable" => "Displacements", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "0")))

    bcs = PeriLab.Solver.Boundary_conditions.boundary_condition(params, test_Data_manager)

    @test length(bcs) == 2
    @test "BC_1" in keys(bcs)
    @test "BC_2" in keys(bcs)

    # params representation
    @test bcs["BC_1"]["Variable"] == "Forces"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]
    @test bcs["BC_2"]["Variable"] == "Displacements"
    @test bcs["BC_2"]["Coordinate"] == "z"
    @test bcs["BC_2"]["Value"] == "0"
    @test bcs["BC_2"]["Node Set"] == [4, 2, 7, 10]
    @test !("BC_3" in keys(bcs))
end
@testset "ut_check_valid_bcs" begin
    test_Data_manager = PeriLab.Data_manager
    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Variable" => "Forces", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Variable" => "not there", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "0")))

    test_Data_manager.set_nset("Nset_1", [1, 2, 3])
    test_Data_manager.set_nset("Nset_2", [3, 4, 7, 10])
    test_Data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 3, 3 => 4, 4 => 2, 5 => 5, 6 => 6, 7 => 7, 8 => 8, 9 => 9, 10 => 10))

    bcs = PeriLab.Solver.Boundary_conditions.boundary_condition(params, test_Data_manager)
    @test isnothing(PeriLab.Solver.Boundary_conditions.check_valid_bcs(bcs, test_Data_manager))
end
@testset "ut_init_BCs" begin

    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.set_num_controller(10)

    test_Data_manager.create_constant_node_field("Coordinates", Float64, 3)
    test_Data_manager.create_constant_node_field("Forces", Float64, 3)
    test_Data_manager.create_node_field("Displacements", Float64, 3)
    test_Data_manager.set_dof(2)

    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Variable" => "Forces", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Variable" => "Displacements", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "5")))

    bcs = PeriLab.Solver.Boundary_conditions.init_BCs(params, test_Data_manager)
    @test length(bcs) == 1
    # clean up params representation
    @test "BC_1" in keys(bcs)
    @test ("BC_2" in keys(bcs)) == false
    @test bcs["BC_1"]["Variable"] == "Forces"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]

    test_Data_manager.set_dof(3)
    bcs = PeriLab.Solver.Boundary_conditions.init_BCs(params, test_Data_manager)
    @test length(bcs) == 2
    @test "BC_1" in keys(bcs)
    @test "BC_2" in keys(bcs)
    @test bcs["BC_1"]["Variable"] == "Forces"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]
    @test bcs["BC_2"]["Variable"] == "DisplacementsNP1"
    @test bcs["BC_2"]["Coordinate"] == "z"
    @test bcs["BC_2"]["Value"] == "5"
    @test bcs["BC_2"]["Node Set"] == [4, 2, 7, 10]

end

@testset "ut_apply_bc" begin

    test_Data_manager = PeriLab.Data_manager

    test_Data_manager.set_dof(3)
    test_Data_manager.create_constant_node_field("Coordinates", Float64, 3)
    test_Data_manager.create_constant_node_field("Forces", Float64, 3)
    test_Data_manager.create_node_field("Displacements", Float64, 3)

    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Variable" => "Forces", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Variable" => "Displacements", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "5")))

    force = test_Data_manager.get_field("Forces")
    disp = test_Data_manager.get_field("Displacements", "NP1")
    @test sum(force) == 0
    @test sum(disp) == 0
    bcs = PeriLab.Solver.Boundary_conditions.init_BCs(params, test_Data_manager)
    PeriLab.Solver.Boundary_conditions.apply_bc_dirichlet(bcs, test_Data_manager, 0.0)
    force = test_Data_manager.get_field("Forces")
    disp = test_Data_manager.get_field("Displacements", "NP1")
    @test sum(force) == 0
    @test sum(disp) == 20
    @test disp == [0 0 0; 0 0 5; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5]

    PeriLab.Solver.Boundary_conditions.apply_bc_dirichlet(bcs, test_Data_manager, 0.2)
    force = test_Data_manager.get_field("Forces")
    disp = test_Data_manager.get_field("Displacements", "NP1")
    @test sum(force) == 12
    @test force == [4 0 0; 0 0 0; 4 0 0; 4 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]
    @test sum(disp) == 20
    @test disp == [0 0 0; 0 0 5; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5]
    # test if global nodes are not at the core
    bcs["BC_1"]["Node Set"] = []
    bcs["BC_2"]["Node Set"] = []
    force .= 0
    disp .= 0
    PeriLab.Solver.Boundary_conditions.apply_bc_dirichlet(bcs, test_Data_manager, 0.2)
    @test sum(force) == 0
    @test sum(disp) == 0

    params = Dict("Boundary Conditions" => Dict("BC_2" => Dict("Variable" => "Displacements", "Node Set" => "Nset_2", "Coordinate" => "u", "Value" => "5")))
    bcs = PeriLab.Solver.Boundary_conditions.init_BCs(params, test_Data_manager)
    @test isnothing(PeriLab.Solver.Boundary_conditions.apply_bc_dirichlet(bcs, test_Data_manager, 0.2))
end