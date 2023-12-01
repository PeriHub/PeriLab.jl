# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Documenter, PeriLab, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
    plugins=[bib],
    modules=[PeriLab],
    authors="Christian Willberg <christian.willberg@dlr.de> and Jan-Timo Hesse <jan-timo.hesse@dlr.de>",
    doctest=true,
    checkdocs=:none, # :all, :exports, :none
    sitename="PeriLab",
    repo=Documenter.Remotes.GitLab("dlr-PeriHub", "PeriLab"),
    format=Documenter.HTML(
        canonical="https://gitlab.com/dlr-perihub/perilab",
        assets=["assets/favicon.ico"],
        edit_link="main"
    ),
    pages=Any[
        "Introduction"=>"index.md",
        "First Steps with PeriLab"=>"man/basics.md",
        "User Guide"=>Any[
            "Getting Started"=>"man/getting_started.md",
            "Input File"=>"man/input_yaml.md",
            "Mesh and Nodesets"=>"man/mesh_input.md",
        ],
        "API"=>Any[
            # "Types" => "lib/types.md",
            "Functions"=>"lib/functions.md",
            "IO Functions"=>"lib/io_functions.md",
            "Physics Functions"=>"lib/physics_functions.md",
            # "Indexing" => "lib/indexing.md",
            # "Metadata" => "lib/metadata.md",
            # hide("Internals" => "lib/internals.md"),
        ],
        "References"=>"lib/references.md",
        "Dev Log"=>"devLog.md",
    ]
)
