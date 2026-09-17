# Module overview

The following diagram shows how the main modules of `PeriLab.jl` fit together:
the CLI entry, the `run()` orchestrator, the central `Data_Manager`, the
`solver()` time-stepping loop, the pluggable multi-physics `Model` factory, the
`FEM` coupling, and the shared support/compute/MPI foundations.

```@raw html
<style>
.archify-view {
    display: block;
    width: 100%;
    border: 0;
    max-width: 1280px;
}
</style>
<iframe src="assets/PeriLab-overview.html" class="archify-view" height="720" frameborder="0"></iframe>
```

The module hierarchy below lists the submodules in more detail.

- PeriLab
    - Helpers
    - Geometry
    - Data_Manager
    - Logging_Module
    - MPI_Communication
    - Parameter_Handling
    - IO
    - Solver_Manager
        - Material_Basis
        - FEM
            - FEM_Basis
            - Coupling
                - Arlequin_Coupling
        - Model
            - Pre_Calculation
                - Axissymmetric
                - Bond_Deformation
                - Deformation_Gradient
                - Pre_Bond_Associated_Correspondence
                - Shape_Tensor
                - …
            - Surface_Correction
            - Contact
                - Contact_Search
                - Penalty_Model
                - …
            - Additive
                - Damage_Based
                - …
            - Degradation
                - Thermal_Decomposition
                - …
            - Damage
                - Critical_Stretch
                - Critical_Energy
                - Critical_Energy_Aniso
                - ...
            - Material
                - Ordinary
                - OneD_Bond_Based_Elastic
                - Bondbased_Elastic
                - Unified_Bondbased_Elastic
                - Correspondence
                    - Global_Zero_Energy_Control
                    - Bond_Associated_Correspondence
                    - Correspondence_Elastic
                    - Correspondence_Plastic
                    - Correspondence_UMAT
                    - Correspondence_VUMAT
                    - …
                - PD_Solid_Elastic
                - PD_Solid_Plastic
                - …
            - Thermal
                - Heat_Transfer
                - HETVAL
                - Thermal_Expansion
                - Thermal_Flow
                - ...
        - Boundary_Conditions
        - Verlet_Solver
        - Static_Solver
        - Influence_Function
