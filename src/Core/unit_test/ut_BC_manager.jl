include("../BC_manager.jl")
include("../../Support/data_manager.jl")
using Test
import .Data_manager
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

@testset "ut_boundary_condition" begin
    testDatamanager = Data_manager
    testDatamanager.set_dof() = 2
    params = Dict()
    bcs = Boundary_conditions.boundary_condition(params, testDatamanager)
    @test length(bcs) == 0
    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "0")))

    testDatamanager.set_nsets("Nset_1", [1, 2, 3])
    testDatamanager.set_nsets("Nset_2", [3, 4, 7, 10])
    testDatamanager.set_glob_to_loc([1, 3, 4, 2, 5, 6, 7, 8, 9, 10])

    bcs = Boundary_conditions.boundary_condition(params, testDatamanager)
    @test length(bcs) == 2
    @test "BC_1" in keys(bcs)
    @test "BC_2" in keys(bcs)
    # params representation
    @test bcs["BC_1"]["Type"] == "Force"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]
    @test bcs["BC_2"]["Type"] == "Displacement"
    @test bcs["BC_2"]["Coordinate"] == "z"
    @test bcs["BC_2"]["Value"] == "0"
    @test bcs["BC_2"]["Node Set"] == [4, 2, 7, 10]
end

@testset "ut_init_BCs" begin

    testDatamanager = Data_manager
    testDatamanager.set_nnodes(10)

    testDatamanager.create_constant_node_field("Coordinates", Float32, 3)
    testDatamanager.create_constant_node_field("Force", Float32, 3)
    testDatamanager.create_node_field("Displacements", Float32, 3)
    testDatamanager.set_dof(2)

    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "5")))

    bcs = Boundary_conditions.init_BCs(params, testDatamanager)
    @test length(bcs) == 1
    # clean up params representation
    @test "BC_1" in keys(bcs)
    @test ("BC_2" in keys(bcs)) == false
    @test bcs["BC_1"]["Type"] == "Force"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]

    testDatamanager.set_dof(3)
    bcs = Boundary_conditions.init_BCs(params, testDatamanager)
    @test length(bcs) == 2
    @test "BC_1" in keys(bcs)
    @test "BC_2" in keys(bcs)
    @test bcs["BC_1"]["Type"] == "Force"
    @test bcs["BC_1"]["Coordinate"] == "x"
    @test bcs["BC_1"]["Value"] == "20*t"
    @test bcs["BC_1"]["Node Set"] == [1, 3, 4]
    @test bcs["BC_2"]["Type"] == "DisplacementsNP1"
    @test bcs["BC_2"]["Coordinate"] == "z"
    @test bcs["BC_2"]["Value"] == "5"
    @test bcs["BC_2"]["Node Set"] == [4, 2, 7, 10]

end

@testset "ut_apply_bc" begin

    testDatamanager = Data_manager

    testDatamanager.reset_filter()
    testDatamanager.set_dof(3)
    testDatamanager.create_constant_node_field("Coordinates", Float32, 3)
    testDatamanager.create_constant_node_field("Force", Float32, 3)
    testDatamanager.create_node_field("Displacements", Float32, 3)
    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "z", "Value" => "5")))

    force = testDatamanager.get_field("Force")
    disp = testDatamanager.get_field("Displacements", "NP1")
    @test sum(force) == 0
    @test sum(disp) == 0
    bcs = Boundary_conditions.init_BCs(params, testDatamanager)
    Boundary_conditions.apply_bc(bcs, testDatamanager, 0)
    force = testDatamanager.get_field("Force")
    disp = testDatamanager.get_field("Displacements", "NP1")
    @test sum(force) == 0
    @test sum(disp) == 20
    @test disp == [0 0 0; 0 0 5; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5]

    Boundary_conditions.apply_bc(bcs, testDatamanager, 0.2)
    force = testDatamanager.get_field("Force")
    disp = testDatamanager.get_field("Displacements", "NP1")
    @test sum(force) == 12
    @test force == [4 0 0; 0 0 0; 4 0 0; 4 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]
    @test sum(disp) == 20
    @test disp == [0 0 0; 0 0 5; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5; 0 0 0; 0 0 0; 0 0 5]
end