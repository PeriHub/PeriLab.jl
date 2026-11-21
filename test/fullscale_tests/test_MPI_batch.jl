# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# MPI runs of the fullscale decks, batched: one mpiexec per rank count, so PeriLab is
# compiled once per batch instead of once per deck (40-90 s per launch for decks that run
# in seconds). Each deck writes into <folder>/mpi_<ranks> and is compared with its serial
# reference like run_perilab does.
#   (folder, deck, ranks, compare_csv)
const MPI_DECKS = [("test_PD_solid_elastic_3D", "strain_xy", 3, false),
    ("test_global_variables", "global_variables", 2, true),
    ("test_Bond_Based_Elastic", "test_bond_based_elastic_4_points", 2, false),
    ("test_Bond_Based_Elastic", "test_bond_based_elastic_3D", 3, false),
    ("test_Penalty_Contact", "penalty_contact_options", 2, false),
    ("test_Penalty_Contact", "penalty_contact_options", 3, false),
    # damage only in the right quarter of the strip: on 2 ranks one rank holds no damage
    # block, but all ranks must still exchange the damage model's fields
    ("test_MPI_partial_damage", "partial_damage", 2, false),
    # stops at the first damage, which only one rank has
    ("test_MPI_partial_damage", "partial_damage_max", 2, false),
    ("test_bond_value_export", "bond_value_export", 2, false),
    ("test_bond_value_export", "bond_value_export_block_wise", 2, false)]

for ranks in unique(d[3] for d in MPI_DECKS)
    decks = [(abspath("fullscale_tests", f), d) for (f, d, r, _) in MPI_DECKS if r == ranks]
    # a failing deck must not stop the batch: report it and exit non-zero at the end
    script = """
    import PeriLab
    failed = false
    for (dir, deck) in $(repr(decks))
        try
            cd(() -> PeriLab.run(deck * ".yaml"; silent = true, output_dir = "mpi_$ranks"), dir)
        catch e
            global failed = true
            @error "MPI batch: \$deck failed" exception = (e, catch_backtrace())
        end
    end
    exit(failed ? 1 : 0)
    """
    @testset "$ranks ranks" begin
        @test success(run(ignorestatus(`$(mpiexec()) -n $ranks $(Base.julia_cmd()) -e $script`)))
        for (folder, deck, r, csv) in MPI_DECKS
            r == ranks || continue
            @testset "$deck" begin
                cd("fullscale_tests/" * folder) do
                    compare_results(deck, folder, ranks; compare_csv = csv,
                                    dir = "mpi_$ranks")
                end
            end
        end
        # keep an output folder only for failed comparisons (their .e and exodiff log)
        for folder in unique(f for (f, _, r, _) in MPI_DECKS if r == ranks)
            dir = joinpath("fullscale_tests", folder, "mpi_$ranks")
            any(f -> endswith(f, ".e") || endswith(f, "_exodiff.log"), readdir(dir)) ||
                rm(dir; recursive = true)
        end
    end
end
