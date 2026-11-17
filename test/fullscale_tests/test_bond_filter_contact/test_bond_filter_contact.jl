# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Two correspondence blocks across a bond filter with `Allow Contact`: the faces must
# transmit compression, slide freely && open freely.

const NX, NY = 12, 16 # nodes; blocks y <= 7 && y >= 8, filter plane at y = 7.5

function write_contact_deck(name, ux, uy; filter = true,
                            material = "Correspondence Elastic")
    bond_filter = filter ? """
    Bond Filters:
      crack:
        Type: Rectangular_Plane
        Normal X: 0.0
        Normal Y: 1.0
        Normal Z: 0.0
        Lower Left Corner X: -1.0
        Lower Left Corner Y: 7.5
        Lower Left Corner Z: -1.0
        Bottom Unit Vector X: 1.0
        Bottom Unit Vector Y: 0.0
        Bottom Unit Vector Z: 0.0
        Bottom Length: $(NX + 2.0)
        Side Length: 2.0
        Allow Contact: true
""" : ""
    write(name * ".yaml",
          """
          PeriLab:
            Discretization:
              Type: Text File
              Input Mesh File: contact_mesh.txt
              Node Sets:
                bottom: "y < 2.5"
                top: "y > 12.5"
                lower_face: "y > 6.5 && y < 7.5"
                upper_face: "y > 7.5 && y < 8.5"
          $bond_filter
            Models:
              Pre Calculation Global:
                Deformed Bond Geometry: true
                Deformation Gradient: true
                Shape Tensor: true
              Material Models:
                Mat:
                  Material Model: $material
                  Symmetry: isotropic plane stress
                  Young's Modulus: 1000.0
                  Poisson's Ratio: 0.3
          $(material == "Correspondence Elastic" ? "        Zero Energy Control: Global" : "")
              Damage Models:
                Damage:
                  Damage Model: Critical Stretch
                  Critical Value: 10.0
            Blocks:
              block_1:
                Block ID: 1
                Material Model: Mat
                Damage Model: Damage
                Density: 1.0e-6
                Horizon: 3.015
              block_2:
                Block ID: 2
                Material Model: Mat
                Damage Model: Damage
                Density: 1.0e-6
                Horizon: 3.015
            Boundary Conditions:
              BC_1:
                Variable: Displacements
                Node Set: bottom
                Coordinate: x
                Value: 0.0
                Type: Dirichlet
              BC_2:
                Variable: Displacements
                Node Set: bottom
                Coordinate: y
                Value: 0.0
                Type: Dirichlet
              BC_3:
                Variable: Displacements
                Node Set: top
                Coordinate: x
                Value: $ux*t
                Type: Dirichlet
              BC_4:
                Variable: Displacements
                Node Set: top
                Coordinate: y
                Value: $uy*t
                Type: Dirichlet
            Compute Class Parameters:
              F_top:
                Compute Class: Node_Set_Data
                Calculation Type: Sum
                Node Set: top
                Variable: Forces
              F_bottom:
                Compute Class: Node_Set_Data
                Calculation Type: Sum
                Node Set: bottom
                Variable: Forces
              u_lower:
                Compute Class: Node_Set_Data
                Calculation Type: Maximum
                Node Set: lower_face
                Variable: Displacements
              u_upper:
                Compute Class: Node_Set_Data
                Calculation Type: Minimum
                Node Set: upper_face
                Variable: Displacements
              damage_face:
                Compute Class: Node_Set_Data
                Calculation Type: Minimum
                Node Set: lower_face
                Variable: Damage
            Solver:
              Initial Time: 0.0
              Final Time: 1.0
              Number of Steps: 4
              Damage Models: true
              Static:
                Maximum number of iterations: 100
                Residual tolerance: 1e-8
                Solution tolerance: 1e-8
                Residual scaling: 1000
                Linear Start Value: 0 0 0 0
            Outputs:
              Output1:
                Output File Type: CSV
                Output Filename: $name
                Output Frequency: 1
                Output Variables:
                  F_top: true
                  F_bottom: true
                  u_lower: true
                  u_upper: true
                  damage_face: true
          """)
end

# last row of the CSV output as Dict column => value
function last_csv_row(name)
    lines = readlines(name * ".csv")
    header = strip.(split(lines[1], ","))
    values = parse.(Float64, split(lines[end], ","))
    return Dict(zip(header, values))
end

function run_contact(name, ux, uy; kwargs...)
    write_contact_deck(name, ux, uy; kwargs...)
    PeriLab.run(name * ".yaml"; silent = true)
    result = last_csv_row(name)
    rm(name * ".yaml")
    rm(name * ".csv")
    return result
end

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    open("contact_mesh.txt", "w") do io
        println(io, "header: x y block_id volume")
        for y in 0:(NY - 1), x in 0:(NX - 1)
            println(io, "$x.0 $y.0 $(y < 8 ? 1 : 2) 1.0")
        end
    end
    d = 0.05
    # approach of the face rows (gap < 0) and the load in an intact model
    glued = run_contact("glued", 0.0, -d; filter = false)
    glued_gap = glued["u_uppery"] - glued["u_lowery"]

    for material in ("Correspondence Elastic", "Bond-based Elastic")
        compression = run_contact("compression", 0.0, -d; material)
        # the load is transmitted
        @test compression["F_topy"] > 0
        @test isapprox(compression["F_topy"], -compression["F_bottomy"]; rtol = 1e-4)
        # the faces only overlap by about the elastic approach of the intact model
        @test compression["u_uppery"] - compression["u_lowery"] > 2 * glued_gap
        # contact bonds count as broken in the damage index
        @test compression["damage_face"] > 0.3
        material == "Correspondence Elastic" &&
            @test isapprox(compression["F_topy"], glued["F_topy"]; rtol = 0.3)
    end

    # frictionless: no shear force across the plane
    glued_shear = run_contact("glued_shear", d, 0.0; filter = false)
    shear = run_contact("shear", d, 0.0)
    @test abs(shear["F_topx"]) < 1e-4 * abs(glued_shear["F_topx"])

    # opening: no force across the plane
    opening = run_contact("opening", 0.0, d)
    @test abs(opening["F_topy"]) < 1e-4 * glued["F_topy"]
    rm("contact_mesh.txt")
end
