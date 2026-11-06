# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const Gcode_Mesh = PeriLab.IO.Mesh_Import.Gcode_Mesh

# mesh (sampling 1 mm, one node per cross section) of the given g-code lines
function gcode_mesh(lines::Vector{String})
    file = tempname() * ".gcode"
    write(file, join(lines, "\n"))
    pd_mesh = Gcode_Mesh.new_pd_mesh(1.0, 0.4, 0.4, 1, 1, nothing)
    commands = Dict{String,Any}("Start" => nothing, "Stop" => nothing, "End" => nothing)
    Gcode_Mesh.write_mesh(file, commands, true, pd_mesh)
    rm(file)
    return pd_mesh["mesh_df"]
end

radius(df, cx, cy) = hypot.(df.x .- cx, df.y .- cy)

# feed 600 mm/min = 10 mm/s; the points lie on chords of the arc (sagitta < 0.02 mm)
@testset "ut_gcode_arc_ij" begin
    # G3: counter-clockwise quarter circle around the origin, (10,0) -> (0,10)
    df = gcode_mesh(["G0 X10 Y0 Z0.2 F600", "G3 X0 Y10 I-10 J0 E1 F600"])
    @test all(abs.(radius(df, 0, 0) .- 10) .< 0.02)
    @test all(df.x .> -1e-9) && all(df.y .> -1e-9)
    @test size(df, 1) == 16 # 15.7 mm arc
    @test issorted(df.Activation_Time)
    @test isapprox(maximum(df.Activation_Time), 2.5; atol = 0.1) # 1 s travel + 1.57 s arc
    # G2: clockwise, (10,0) -> (0,-10) stays in the fourth quadrant
    df = gcode_mesh(["G0 X10 Y0 Z0.2 F600", "G2 X0 Y-10 I-10 J0 E1"])
    @test all(df.x .> -1e-9) && all(df.y .< 1e-9)
    # full circle: end point = start point
    df = gcode_mesh(["G0 X10 Y0 Z0.2 F600", "G3 X10 Y0 I-10 J0 E1"])
    @test size(df, 1) == 63 # 62.8 mm
    @test minimum(df.x) < -9.9
end

@testset "ut_gcode_arc_r" begin
    # positive R: the short arc (<= 180°)
    df = gcode_mesh(["G0 X0 Y10 Z0.2 F600", "G2 X10 Y0 R10 E1"])
    @test all(abs.(radius(df, 0, 0) .- 10) .< 0.02)
    @test size(df, 1) == 16
    df = gcode_mesh(["G0 X0 Y10 Z0.2 F600", "G3 X10 Y0 R10 E1"])
    @test all(abs.(radius(df, 10, 10) .- 10) .< 0.02)
    @test size(df, 1) == 16
    # negative R: the long arc
    df = gcode_mesh(["G0 X0 Y10 Z0.2 F600", "G2 X10 Y0 R-10 E1"])
    @test all(abs.(radius(df, 10, 10) .- 10) .< 0.02)
    @test size(df, 1) == 47 # 47.1 mm
    # modal: a line of only axis words repeats the last arc
    df = gcode_mesh(["G0 X10 Y0 Z0.2 F600", "G3 X0 Y10 I-10 J0 E1", "X-10 Y0 I0 J-10 E2"])
    @test all(abs.(radius(df, 0, 0) .- 10) .< 0.02)
    @test minimum(df.x) < -9.9
end

@testset "ut_gcode_dwell_and_rotation" begin
    # 10 mm at 10 mm/s, dwell 2 s (S in seconds), 10 mm, dwell 0.5 s (P in ms), 10 mm
    df = gcode_mesh(["G0 X0 Y0 Z0.2 F600", "G1 X10 E1", "G4 S2", "G1 X20 E2", "G4 P500",
                        "G1 X30 E3"])
    t(x) = df.Activation_Time[argmin(abs.(df.x .- x))]
    @test isapprox(t(5), 0.5; atol = 0.1)
    @test isapprox(t(15), 3.5; atol = 0.1)
    @test isapprox(t(25), 5.0; atol = 0.1)
    # G1 with C: rotation about the origin at constant radius, (10,0) -> (0,10)
    df = gcode_mesh(["G0 X10 Y0 Z0.2 F600", "G1 C90 E1"])
    @test all(abs.(radius(df, 0, 0) .- 10) .< 0.02)
    @test size(df, 1) == 16
end
