# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const ID = PeriLab.InputDeck

const UT_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))

# Decks that are known not to parse, with the reason. Fix or document them;
# never add a deck here to hide a validator bug.
const UT_DECK_ALLOWLIST = Dict("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml" => "uses a non-existent \"External\" solver; not run by any test")

# licensed models (`ut_license_error`, test/helper.jl): without a license a deck
# naming them cannot be validated

function ut_all_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_REPO_ROOT, root))
        for file in files
            endswith(file, ".yaml") && push!(decks, joinpath(dir, file))
        end
    end
    return sort!(decks)
end

@testset "every shipped deck parses in strict mode" begin
    decks = ut_all_decks()
    @test length(decks) > 100
    for file in decks
        relative = relpath(file, UT_REPO_ROOT)
        haskey(UT_DECK_ALLOWLIST, relative) && continue
        raw = PeriLab.IO.read_input(file)
        if !(raw isa AbstractDict && haskey(raw, "PeriLab"))
            continue                        # not an input deck (e.g. a config file)
        end
        _, ctx = ID.read_input(raw["PeriLab"], dirname(file); strict = true)
        errors = filter(e -> e.severity == :error && !ut_license_error(e), ctx.errors)
        if !isempty(errors)
            @error "$relative\n" * PeriLab.ParameterSpec.format_errors(errors)
        end
        @test isempty(errors)
    end
end
