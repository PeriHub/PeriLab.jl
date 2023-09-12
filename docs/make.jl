using Documenter, PeriLab, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
    bib,
    modules=[PeriLab],
    authors="Christian Willberg <christian.willberg@dlr.de> and Jan-Timo Hesse <jan-timo.hesse@dlr.de>",
    checkdocs=:exports,
    sitename="PeriLab",
    repo="https://gitlab.dlr.de/fa_sw/peridynamik/perilab.git/blob/{commit}{path}#{line}",
    pages=[
        "index.md",
        "Dev Log" => "devLog.md",
        "Compute" => [
            "compute_forces" => "Compute/compute_forces.md"
        ],
        "Core" => [
            "Data_model" => [
                "data" => "Core/Data_model/data.md"
            ],
            "Solver" => [
                "Solver_control" => "Core/Solver/Solver_control.md",
                "Verlet" => "Core/Solver/Verlet.md"
            ],
            "BC_manager" => "Core/BC_manager.md"
        ],
        "IO" => [
            "exodus_export" => "IO/exodus_export.md",
            "IO" => "IO/IO.md",
            "mesh_data" => "IO/mesh_data.md",
            "read_inputdeck" => "IO/read_inputdeck.md",
        ],
        "MPI_communication" => [
            "MPI_communication" => "MPI_communication/MPI_communication.md"
            "MPI_init" => "MPI_communication/MPI_init.md"
        ],
        "Physics" => [
            "Physics" => "Physics/Physics_Factory.md",
            "Additive" => [
                "Additve_Factory" => "Physics/Additive/Additive_Factory.md"
            ],
            "Damage" => [
                "Damage_Factory" => "Physics/Damage/Damage_Factory.md"
            ],
            "Material" => [
                "material_basis" => "Physics/Material/material_basis.md",
                "Material_Factory" => "Physics/Material/Material_Factory.md",
                "BondBased" => [
                    "BondBased" => "Physics/Material/BondBased/Bondbased.md",
                    "Elastic" => "Physics/Material/BondBased/Elastic.md"
                ],
                "Correspondence" => [
                    "Correspondence" => "Physics/Material/Correspondence/Correspondence.md"
                ],
                "Ordinary" => [
                    "Ordinary" => "Physics/Material/Ordinary/Ordinary.md",
                    "PD_Solid_Elastic" => "Physics/Material/Ordinary/PD_Solid_Elastic.md"
                ]
            ],
            "Thermal" => [
                "Thermal_Factory" => "Physics/Thermal/Thermal_Factory.md",
                "BondBased" => [
                    "Thermal_bond_based" => "Physics/Thermal/BondBased/Thermal_bond_based.md",
                ],
                "Correspondence" => [
                    "Thermal_correspondence" => "Physics/Thermal/Correspondence/Thermal_correspondence.md"
                ],
            ],
        ],
        "Support" => [
            "data_manager" => "Support/data_manager.md",
            "geometry" => "Support/geometry.md",
            "helpers" => "Support/helpers.md",
            "tools" => "Support/tools.md",
            "Parameters" => [
                "parameter_handling" => "Support/Parameters/parameter_handling.md",
                "parameter_handling_bc" => "Support/Parameters/parameter_handling_bc.md",
                "parameter_handling_blocks" => "Support/Parameters/parameter_handling_blocks.md",
                "parameter_handling_mesh" => "Support/Parameters/parameter_handling_mesh.md",
                "parameter_handling_output" => "Support/Parameters/parameter_handling_output.md",
                "parameter_handling_physics" => "Support/Parameters/parameter_handling_physics.md",
                "parameter_handling_solver" => "Support/Parameters/parameter_handling_solver.md",
            ]
        ]
    ]
)

# makedocs(sitename="PeriLab")
