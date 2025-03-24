# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Documenter, PeriLab, DocumenterCitations, DocumenterMermaid

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(plugins = [bib],
         modules = [PeriLab],
         authors = "Christian Willberg <christian.willberg@dlr.de> and Jan-Timo Hesse <jan-timo.hesse@dlr.de>",
         doctest = true,
         checkdocs = :none, # :all, :exports, :none
         sitename = "PeriLab",
         repo = Documenter.Remotes.GitHub("PeriHub", "PeriLab.jl"),
         format = Documenter.HTML(canonical = "https://github.com/PeriHub/PeriLab.jl",
                                  assets = ["assets/favicon.ico"],
                                  edit_link = "main"),
         pages = Any["Introduction" => "index.md",
                     "First Steps with PeriLab" => "man/basics.md",
                     "User Guide" => Any["Getting Started" => "man/getting_started.md",
                                         "Input File" => "man/input_yaml.md",
                                         "Mesh and Nodesets" => "man/mesh_input.md",
                                         "Bond-Filter" => "man/bond_filter.md",
                                         "Models" => Any["Overview" => "man/models/overview.md",
                                                         "Material Models" => "man/models/materials.md",
                                                         "Damage Models" => "man/models/damage.md",
                                                         "Thermal Models" => "man/models/thermal.md",
                                                         "Additive Models" => "man/models/additive.md"],
                                         "Output" => "man/output.md",
                                         "Solver" => "man/solver/solver.md"],
                     "Theory" => Any["Basics" => "theory/theory.md",
                                     "Peridynamic Basics" => Any["Bond-based" => "theory/theory_bondbased.md",
                                                                 "Ordinary state-based" => "theory/theory_ordinary.md",
                                                                 "Non-Ordinary state-based" => "theory/theory_correspondence.md"],
                                     "Multi-Physics" => Any["Thermal Models" => "theory/theory_thermal.md"],
                                     "Finite Element Method" => Any["FEM" => "theory/theory_FEM.md",
                                                                    "Finite Element - Peridynamics coupling" => "theory/theory_FEM_PD_coupling.md"],
                                     "Modelling strategies" => Any["Surface Correction" => "theory/theory_surface_correction.md"]],
                     "Lecture" => Any["Non-local structural mechanics and peridynamics" => "lecture/lecture.md",
                                      "Seminar_1" => "lecture/seminar_1.md"],
                     "FAQ" => "lib/faq.md",
                     "Glossar" => "lib/glossar.md",
                     "References" => "lib/references.md",
                     "Useful Links" => "lib/links.md",
                     # "Dev Log"=>"devLog.md",
                     "Development" => Any["Guide" => "man/dev/developement_guide.md",
                                          "Datamanager" => "man/dev/datamanager.md",
                                          "Module integration" => "man/dev/module_integration.md",
                                          "Software Testing" => "man/dev/testing.md"],
                     "API" => Any["Model Factory" => "lib/model_factory_functions.md",
                                  "Data Manager" => "lib/data_manager_functions.md",
                                  "Solver" => "lib/solver_functions.md",
                                  "IO" => "lib/io_functions.md",
                                  "Geometry" => "lib/geometry_functions.md",
                                  "MPI" => "lib/mpi_functions.md",
                                  "Helper" => "lib/helper_functions.md",
                                  "Logging" => "lib/logging_functions.md"
                                  # hide("Internals" => "lib/internals.md"),
                                  ]])
deploydocs(repo = "github.com/PeriHub/PeriLab.jl.git")
