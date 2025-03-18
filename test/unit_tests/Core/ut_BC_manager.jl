# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


using Test

#include("../../../src/PeriLab.jl")
#using .PeriLab
@testset "ut_clean_up" begin
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("") == ""
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("-") == " .- "
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("+") == " .+ "
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("*") == " .* "
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("/") == " ./ "
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("5e-8") == "5e-8"
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("5.0e-8") == "5.0e-8"
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("5.0e-8-5") == "5.0e-8 .- 5"
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("5.0e-8/5e+5*t") ==
          "5.0e-8 ./ 5e+5 .* t"
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("5-8.0/5e+5") ==
          "5 .- 8.0 ./ 5e+5"
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("2") == "2"
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("2.") == "2."
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("2.0") == "2.0"
    @test PeriLab.Solver_control.Boundary_conditions.clean_up("2.0.-1") == "2.0 .- 1"
end

@testset "ut_evaluation" begin
    unit = ones(Float64, 3)
    dof = 2
    time::Float64 = 2
    step_time::Float64 = 2
    coor = zeros(3, 3)
    bc = Int64(10)
    values = ones(Float64, 3)
    PeriLab.Solver_control.Boundary_conditions.eval_bc!(
        values,
        bc,
        coor,
        time,
        step_time,
        dof,
        false,
    )
    @test (10 * unit) == values
    values = ones(Float64, 3)
    bc = Float64(10)
    PeriLab.Solver_control.Boundary_conditions.eval_bc!(
        values,
        bc,
        coor,
        time,
        step_time,
        dof,
        false,
    )
    @test (10 * unit) == values
    values = ones(Float64, 3)
    bc = Float64(10)
    PeriLab.Solver_control.Boundary_conditions.eval_bc!(
        values,
        bc,
        coor,
        time,
        step_time,
        dof,
        false,
    )
    @test (10 * unit) == values
    values = ones(Float64, 3)
    bc = string(10)
    PeriLab.Solver_control.Boundary_conditions.eval_bc!(
        values,
        bc,
        coor,
        time,
        step_time,
        dof,
        false,
    )
    @test (10 * unit) == values
    values = ones(Float64, 3)
    bc = "x"
    for i = 1:4
        coor[3, 1] = i * i - 2
        PeriLab.Solver_control.Boundary_conditions.eval_bc!(
            values,
            bc,
            coor,
            time,
            step_time,
            dof,
            false,
        )
        @test (coor[:, 1]) == values
        values = ones(Float64, 3)
    end
    dof = 3
    bc = "t"
    PeriLab.Solver_control.Boundary_conditions.eval_bc!(
        values,
        bc,
        coor,
        time,
        step_time,
        dof,
        false,
    )
    @test (time * unit) == values
    values = ones(Float64, 3)
    bc = "t*x"
    PeriLab.Solver_control.Boundary_conditions.eval_bc!(
        values,
        bc,
        coor,
        time,
        step_time,
        dof,
        false,
    )
    @test (time .* coor[:, 1]) == values
    values = ones(Float64, 3)
    for t = 0:4
        bc = "if t>2 0 else 20 end"
        if t > 2
            PeriLab.Solver_control.Boundary_conditions.eval_bc!(
                values,
                bc,
                coor,
                Float64(t),
                Float64(t),
                dof,
                false,
            )
            @test (0.0 * unit) == values
            values = ones(Float64, 3)
        else
            PeriLab.Solver_control.Boundary_conditions.eval_bc!(
                values,
                bc,
                coor,
                Float64(t),
                Float64(t),
                dof,
                false,
            )
            @test (20.0 * unit) == values
            values = ones(Float64, 3)
        end
    end
    for t = 0:2
        bc = "100.0"
        if t == 0
            PeriLab.Solver_control.Boundary_conditions.eval_bc!(
                values,
                bc,
                coor,
                Float64(t),
                Float64(t),
                dof,
                false,
            )
            @test (100.0 * unit) == values
            values = ones(Float64, 3)
            PeriLab.Solver_control.Boundary_conditions.eval_bc!(
                values,
                bc,
                coor,
                Float64(t),
                Float64(t),
                dof,
                true,
            )
            @test (100.0 * unit) == values
            values = ones(Float64, 3)
        else
            PeriLab.Solver_control.Boundary_conditions.eval_bc!(
                values,
                bc,
                coor,
                Float64(t),
                Float64(t),
                dof,
                true,
            )
            @test unit == values
            values = ones(Float64, 3)
        end
    end
    temp_values = values
    PeriLab.Solver_control.Boundary_conditions.eval_bc!(
        values,
        bc,
        Matrix{Float64}(undef, 0, 0),
        time,
        step_time,
        dof,
        false,
    )
    @test values == temp_values
end

@testset "ut_boundary_condition" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof() = 2
    params = Dict()
    bcs = PeriLab.Solver_control.Boundary_conditions.boundary_condition(
        params,
        test_data_manager,
    )
    @test length(bcs) == 0
    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Forces",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
            "BC_2" => Dict(
                "Variable" => "Displacements",
                "Node Set" => "Nset_2",
                "Coordinate" => "z",
                "Value" => "0",
            ),
            "BC_3" => Dict(
                "Variable" => "Displacements",
                "Node Set" => "Nset_3",
                "Coordinate" => "z",
                "Value" => "0",
            ),
        ),
    )

    test_data_manager.set_nset("Nset_1", [1, 2, 3])
    test_data_manager.set_nset("Nset_2", [3, 4, 7, 10])
    test_data_manager.set_glob_to_loc(
        Dict(
            1 => 1,
            2 => 3,
            3 => 4,
            4 => 2,
            5 => 5,
            6 => 6,
            7 => 7,
            8 => 8,
            9 => 9,
            10 => 10,
        ),
    )

    bcs = PeriLab.Solver_control.Boundary_conditions.boundary_condition(
        params,
        test_data_manager,
    )
    @test isnothing(bcs)
    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Forces",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
            "BC_2" => Dict(
                "Variable" => "Displacements",
                "Node Set" => "Nset_2",
                "Coordinate" => "z",
                "Value" => "0",
            ),
        ),
    )

    bcs = PeriLab.Solver_control.Boundary_conditions.boundary_condition(
        params,
        test_data_manager,
    )

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
    test_data_manager = PeriLab.Data_manager
    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Forces",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
            "BC_2" => Dict(
                "Variable" => "Displacements",
                "Node Set" => "Nset_2",
                "Coordinate" => "z",
                "Value" => "0",
            ),
        ),
    )

    test_data_manager.set_nset("Nset_1", [1, 2, 3])
    test_data_manager.set_nset("Nset_2", [3, 4, 7, 10])
    test_data_manager.set_glob_to_loc(
        Dict(
            1 => 1,
            2 => 3,
            3 => 4,
            4 => 2,
            5 => 5,
            6 => 6,
            7 => 7,
            8 => 8,
            9 => 9,
            10 => 10,
        ),
    )

    bcs = PeriLab.Solver_control.Boundary_conditions.boundary_condition(
        params,
        test_data_manager,
    )

    # @test bcs = PeriLab.Solver_control.Boundary_conditions.check_valid_bcs(bcs, test_data_manager)


    test_data_manager.set_dof(2)
    @test isnothing(
        PeriLab.Solver_control.Boundary_conditions.check_valid_bcs(bcs, test_data_manager),
    )
    test_data_manager.set_dof(3)


    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Forces",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
            "BC_2" => Dict(
                "Variable" => "not there",
                "Node Set" => "Nset_2",
                "Coordinate" => "z",
                "Value" => "0",
            ),
        ),
    )

    bcs = PeriLab.Solver_control.Boundary_conditions.boundary_condition(
        params,
        test_data_manager,
    )
    @test isnothing(
        PeriLab.Solver_control.Boundary_conditions.check_valid_bcs(bcs, test_data_manager),
    )
end
@testset "ut_init_BCs" begin

    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(10)
    test_data_manager.set_nset("Nset_1", [1, 2, 3])
    test_data_manager.set_nset("Nset_2", [3, 4, 7, 10])
    test_data_manager.set_glob_to_loc(
        Dict(
            1 => 1,
            2 => 3,
            3 => 4,
            4 => 2,
            5 => 5,
            6 => 6,
            7 => 7,
            8 => 8,
            9 => 9,
            10 => 10,
        ),
    )
    test_data_manager.create_constant_node_field("Coordinates", Float64, 3)
    test_data_manager.create_constant_node_field("Forces", Float64, 3)
    test_data_manager.create_node_field("Displacements", Float64, 3)
    test_data_manager.set_dof(3)

    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Forces",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
        ),
    )

    bcs = PeriLab.Solver_control.Boundary_conditions.init_BCs(params, test_data_manager)
    @test length(bcs) == 1
    # clean up params representation
    @test "BC_1" in keys(bcs)
    @test ("BC_2" in keys(bcs)) == false
    @test bcs["BC_1"]["Variable"] == "Forces"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]

    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Forces",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
            "BC_2" => Dict(
                "Variable" => "Displacements",
                "Node Set" => "Nset_2",
                "Coordinate" => "z",
                "Value" => "5",
            ),
        ),
    )
    bcs = PeriLab.Solver_control.Boundary_conditions.init_BCs(params, test_data_manager)
    @test length(bcs) == 2
    @test "BC_1" in keys(bcs)
    @test "BC_2" in keys(bcs)
    @test bcs["BC_1"]["Variable"] == "Forces"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]
    @test bcs["BC_2"]["Variable"] == "Displacements"
    @test bcs["BC_2"]["Time"] == "NP1"
    @test bcs["BC_2"]["Coordinate"] == "z"
    @test bcs["BC_2"]["Value"] == "5"
    @test bcs["BC_2"]["Node Set"] == [4, 2, 7, 10]

end

@testset "ut_apply_bc" begin

    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(10)
    test_data_manager.set_dof(3)
    test_data_manager.set_step(nothing)
    test_data_manager.create_constant_node_field("Coordinates", Float64, 3)
    test_data_manager.create_constant_node_field("Temperature", Float64, 3)
    test_data_manager.create_node_field("Displacements", Float64, 3)
    test_data_manager.create_constant_node_field("Forces", Float64, 3)
    test_data_manager.create_constant_node_field("External Forces", Float64, 3)
    test_data_manager.create_constant_node_field("Force Densities", Float64, 3)
    test_data_manager.create_constant_node_field("External Force Densities", Float64, 3)
    test_data_manager.create_constant_node_field("Density", Float64, 1)
    test_data_manager.set_nset("Nset_1", [1, 2, 3])
    test_data_manager.set_nset("Nset_2", [3, 4, 7, 10])
    test_data_manager.set_glob_to_loc(
        Dict(
            1 => 1,
            2 => 3,
            3 => 4,
            4 => 2,
            5 => 5,
            6 => 6,
            7 => 7,
            8 => 8,
            9 => 9,
            10 => 10,
        ),
    )
    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Temperature",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
            "BC_2" => Dict(
                "Variable" => "Displacements",
                "Node Set" => "Nset_2",
                "Coordinate" => "z",
                "Value" => "5",
            ),
        ),
    )

    temperature = test_data_manager.get_field("Temperature")
    disp = test_data_manager.get_field("Displacements", "NP1")
    @test sum(temperature) == 0
    @test sum(disp) == 0
    bcs = PeriLab.Solver_control.Boundary_conditions.init_BCs(params, test_data_manager)
    PeriLab.Solver_control.Boundary_conditions.apply_bc_dirichlet(
        ["Displacements", "Temperature"],
        bcs,
        test_data_manager,
        0.0,
        0.0,
    )
    temperature = test_data_manager.get_field("Temperature")
    disp = test_data_manager.get_field("Displacements", "NP1")
    @test sum(temperature) == 0
    @test sum(disp) == 20
    @test isapprox(
        disp,
        [0 0 0; 0 0 5; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5],
    )
    PeriLab.Solver_control.Boundary_conditions.apply_bc_dirichlet(
        ["Displacements", "Temperature"],
        bcs,
        test_data_manager,
        0.2,
        0.2,
    )
    temperature = test_data_manager.get_field("Temperature")
    disp = test_data_manager.get_field("Displacements", "NP1")

    @test sum(temperature) == 12
    @test isapprox(
        temperature,
        [4 0 0; 0 0 0; 4 0 0; 4 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0],
    )
    @test sum(disp) == 20
    @test isapprox(
        disp,
        [0 0 0; 0 0 5; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5],
    )
    # test if global nodes are not at the core
    bcs["BC_1"]["Node Set"] = []
    bcs["BC_2"]["Node Set"] = []
    temperature .= 0
    disp .= 0
    PeriLab.Solver_control.Boundary_conditions.apply_bc_dirichlet(
        ["Displacements", "Temperature"],
        bcs,
        test_data_manager,
        0.2,
        0.2,
    )
    @test sum(temperature) == 0
    @test sum(disp) == 0

    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_2" => Dict(
                "Variable" => "Displacements",
                "Node Set" => "Nset_2",
                "Coordinate" => "u",
                "Value" => "5",
            ),
        ),
    )
    bcs = PeriLab.Solver_control.Boundary_conditions.init_BCs(params, test_data_manager)
    @test isnothing(
        PeriLab.Solver_control.Boundary_conditions.apply_bc_dirichlet(
            ["Displacements", "Temperature"],
            bcs,
            test_data_manager,
            0.2,
            0.2,
        ),
    )

    ### apply_bc_dirichlet_force

    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Forces",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20",
            ),
        ),
    )
    bcs = PeriLab.Solver_control.Boundary_conditions.init_BCs(params, test_data_manager)

    PeriLab.Solver_control.Boundary_conditions.apply_bc_dirichlet(
        ["Forces", "Force Densities"],
        bcs,
        test_data_manager,
        0.0,
        0.0,
    )
    force_densities = test_data_manager.get_field("External Forces")
    @test force_densities == [
        20.0 0.0 0.0
        0.0 0.0 0.0
        20.0 0.0 0.0
        20.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
    ]


    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Force Densities",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20",
            ),
        ),
    )
    bcs = PeriLab.Solver_control.Boundary_conditions.init_BCs(params, test_data_manager)

    PeriLab.Solver_control.Boundary_conditions.apply_bc_dirichlet(
        ["Forces", "Force Densities"],
        bcs,
        test_data_manager,
        0.0,
        0.0,
    )
    ext_force_densities = test_data_manager.get_field("External Force Densities")
    @test ext_force_densities == [
        20.0 0.0 0.0
        0.0 0.0 0.0
        20.0 0.0 0.0
        20.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0
    ]


    ### apply_bc_neumann

    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Density",
                "Node Set" => "Nset_1",
                "Value" => "10",
                "Type" => "Neumann",
            ),
        ),
    )
    bcs = PeriLab.Solver_control.Boundary_conditions.init_BCs(params, test_data_manager)

    PeriLab.Solver_control.Boundary_conditions.apply_bc_neumann(
        bcs,
        test_data_manager,
        0.0,
        0.0,
    )
    density = test_data_manager.get_field("Density")
    @test density == [10.0, 0.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    PeriLab.Solver_control.Boundary_conditions.apply_bc_neumann(
        bcs,
        test_data_manager,
        0.0,
        0.0,
    )
    density = test_data_manager.get_field("Density")
    @test density == [20.0, 0.0, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

end
