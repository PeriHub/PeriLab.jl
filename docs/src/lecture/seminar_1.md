<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

## Lecture 1: Bond based truss solution


## Mesh
```plaintext
header: x y block_id volume
0 0 1 1
1 0 1 1
2 0 1 1
3 0 1 1
4 0 1 1
```

## Yaml
```yaml
PeriLab:
    Discretization:
        Node Sets:
        Node Set 1: 1
        Node Set 1: 5
        Type: "Text File"
        Input Mesh File: "mesh.txt"
    Models:
        Material Models:
            Test:
                Material Model: "Bond-based elastic"
                Symmetry: "isotropic plane stress"
                Young's Modulus: 7000
                Poisson's Ratio: 0.3
    Blocks:
        block_1:
            Block ID: 1
            Material Model: "Test"
            Density: 20
            Horizon: 2
    Boundary Conditions:
        BC_1:
            Variable: "Displacement"
            Node Set: "Node Set 1"
            Coordinate: "x"
            Value: "100*t"
            Type: Dirichlet
        BC_2:
            Variable: "Displacement"
            Coordinate: "x"
            Node Set: "Node Set 2"
            Value: "0.1*t"
            Type: Dirichlet
    Solver:
        Material Models: True
        Initial Time: 0.0
        Final Time: 1.0
        Number of Steps: 20
        Static:
            Show solver iteration: true
            Residual tolerance: 1e-7
            Solution tolerance: 1e-8
            Residual scaling: 7000
            m: 550
            Maximum number of iterations: 100
    Outputs:
        Output1:
        Output Filename: "truss"
        Output File Type: Exodus
        Number of Output Steps: 20
        Output Variables:
            Displacements: True
            Number of Neighbors: True
            Forces: True
```
