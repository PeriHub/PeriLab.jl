using Documenter, PeriLab

makedocs(
    modules=[PeriLab],
    authors="Christian Willberg <christian.willberg@dlr.de> and Jan-Timo Hesse <jan-timo.hesse@dlr.de>",
    checkdocs=:exports,
    sitename="PeriLab",
    repo="https://gitlab.dlr.de/fa_sw/peridynamik/perilab.git/blob/{commit}{path}#{line}",
    pages=[
        "index.md",
        "Dev Log" => "devLog.md",
        "Core" => [
            "data_manager" => "data_manager.md"
        ],
        "Physics" => [
            "Material" => [
                "BondBased" => [
                    "elastic" => "elastic.md"
                ]
            ]
        ]
    ]
)

# makedocs(sitename="PeriLab")
