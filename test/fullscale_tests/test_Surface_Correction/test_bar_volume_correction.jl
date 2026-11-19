# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Bar on a regular grid (spacing 1) under a homogeneous uniaxial stress state: every node
# is prescribed, one explicit step. The modulus is the mean PD stress sigma_xx / eps.
function surface_correction_bar(dir, dof, material, correction)
    n = dof == 2 ? (40, 16, 1) : (16, 14, 14)
    horizon = dof == 2 ? 3.015 : 2.015
    nu = material == "bb" ? (dof == 2 ? 1 / 3 : 1 / 4) : 0.3
    eps = 1.0e-3
    cd(dir) do
        open("bar.txt", "w") do io
            println(io,
                    dof == 2 ? "header: x y block_id volume" :
                    "header: x y z block_id volume")
            for k in 1:n[3], j in 1:n[2], i in 1:n[1]
                println(io,
                        join(dof == 2 ? (i - 1, j - 1, 1, 1.0) :
                             (i - 1, j - 1, k - 1, 1, 1.0), " "))
            end
        end
        open("bar_ns.txt", "w") do io
            println(io, "header: global_id")
            foreach(i -> println(io, i), 1:prod(n))
        end
        mat = material == "bb" ?
              ["Material Model: \"Bond-based Elastic\""] :
              ["Material Model: \"Correspondence Elastic\"", "Poisson's Ratio: $nu",
            "Zero Energy Control: \"Global\""]
        push!(mat, "Symmetry: \"isotropic plane stress\"", "Young's Modulus: 1000.0")
        bcs = join(["    BC_$c:\n      Variable: Displacements\n      Node Set: all\n" *
                    "      Coordinate: $c\n      Value: \"$(c == "x" ? eps : -nu * eps)*$c\"\n" *
                    "      Type: Dirichlet\n" for c in ("x", "y", "z")[1:dof]])
        sc = correction ? "  Surface Correction:\n    Type: Volume Correction\n" : ""
        write("bar.yaml",
              """
              PeriLab:
                Discretization:
                  Node Sets:
                    all: bar_ns.txt
                  Type: Text File
                  Input Mesh File: bar.txt
                Models:
                  Material Models:
                    Mat:
                      $(join(mat, "\n        "))
              $(sc)  Blocks:
                  block_1:
                    Block ID: 1
                    Material Model: Mat
                    Density: 1.0
                    Horizon: $horizon
                Boundary Conditions:
              $(bcs)  Solver:
                  Initial Time: 0.0
                  Final Time: 1.0e-6
                  Number of Steps: 1
                  Calculate Cauchy: true
                  Verlet:
                    Fixed dt: 1.0e-6
                Outputs:
                  Output1:
                    Output Filename: bar
                    Output File Type: Exodus
                    Output Frequency: 1
                    Output Variables:
                      Forces: true
              """)
        PeriLab.run("bar.yaml"; silent = true)
    end
    DM = PeriLab.Data_Manager
    x = DM.get_field("Coordinates")
    stress = DM.get_field("Cauchy Stress", "N")[:, 1, 1] ./ eps
    # nodes whose neighbors all have a full horizon: their bonds get the factor 1
    inner = [all(2 * horizon < x[iID, d] < n[d] - 1 - 2 * horizon for d in 1:dof)
             for iID in axes(x, 1)]
    corr = correction ? deepcopy(DM.get_field("Volume Correction")) : nothing
    return (; E_eff = sum(stress) / length(stress),
            E_bulk = sum(stress[inner]) / count(inner), inner, corr,
            force = copy(DM.get_field("Force Densities", "N")), stress)
end

@testset "bar volume correction" begin
    dir = mktempdir()
    for dof in (2, 3)
        plain = surface_correction_bar(dir, dof, "bb", false)
        corrected = surface_correction_bar(dir, dof, "bb", true)
        @test any(plain.inner)
        factors = reduce(vcat, collect.(corrected.corr[findall(plain.inner)]))
        @test all(abs.(factors .- 1) .< 1e-12)
        @test isapprox(plain.E_bulk, 1000; rtol = 0.06)
        @test corrected.E_bulk ≈ plain.E_bulk
        # free surfaces soften the bar; the correction reduces the error of the mean
        # modulus to the bulk one (it overshoots at surfaces parallel to the load: the
        # factor is isotropic, the missing bonds are not)
        @test plain.E_eff < plain.E_bulk
        @test abs(corrected.E_eff - corrected.E_bulk) < abs(plain.E_eff - plain.E_bulk)

        # correspondence: the shape tensor accounts for missing neighbors, not corrected
        plain = surface_correction_bar(dir, dof, "corr", false)
        corrected = surface_correction_bar(dir, dof, "corr", true)
        @test corrected.force == plain.force
        @test corrected.stress == plain.stress
        @test !any(PeriLab.Data_Manager.get_field("Surface Corrected"))
    end
    # scaling the bond forces is a typed loop: no allocation per bond
    sc = PeriLab.Solver_Manager.Model_Factory.Surface_Correction
    PeriLab.Data_Manager.get_field("Surface Corrected") .= true
    scale(nodes) = @allocated sc.compute_surface_volume_correction(nodes)
    scale(1:1)
    @test scale(1:PeriLab.Data_Manager.get_nnodes()) == scale(1:1)
end
