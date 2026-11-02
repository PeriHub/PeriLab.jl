# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# every module an API page lists in `@autodocs Modules = [...]` must exist (docs/make.jl
# runs with CurrentModule = PeriLab and fails on an unknown module)
const UT_DOCS_SRC = normpath(joinpath(@__DIR__, "..", "..", "docs", "src"))

function ut_resolve_module(path::AbstractString)
    m = PeriLab
    for part in split(path, ".")
        isdefined(m, Symbol(part)) || return nothing
        m = getfield(m, Symbol(part))
    end
    return m
end

@testset "docs @autodocs modules exist" begin
    for (dir, _, files) in walkdir(UT_DOCS_SRC), file in files
        endswith(file, ".md") || continue
        for line in eachline(joinpath(dir, file))
            m = match(r"^Modules\s*=\s*\[(.*)\]", strip(line))
            m === nothing && continue
            for name in strip.(split(m.captures[1], ","))
                resolved = ut_resolve_module(name)
                if !(resolved isa Module)
                    @error "$(relpath(joinpath(dir, file), UT_DOCS_SRC)): unknown module $name"
                end
                @test resolved isa Module
            end
        end
    end
end
