# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export Dependent, Constant, Table1D, value, combine

"""
    Dependent

A parameter that is either a constant (`Constant`) or depends on one node
field via a data table (`Table1D`). Read it with `value(d, iID)`.
"""
abstract type Dependent end

struct Constant <: Dependent
    value::Float64
end

"""
    Table1D

Values interpolated (spline) over one node field, read from a data file.
Must be bound to the field array (`bind_table!`) before `value` is called.
"""
struct Table1D <: Dependent
    field_name::String
    x::Vector{Float64}
    y::Vector{Float64}
    spline::Spline1D
    field::Base.RefValue{Vector{Float64}}
    bound::Base.RefValue{Bool}
    warn::Base.RefValue{Bool}
    source::String
end

function Table1D(field_name::AbstractString, x::Vector{Float64}, y::Vector{Float64},
                 source::AbstractString)
    k = min(3, length(x) - 1)
    return Table1D(String(field_name), x, y, Spline1D(x, y; k = k, bc = "nearest"),
                   Ref(Float64[]), Ref(false), Ref(true), String(source))
end

@inline value(c::Constant, ::Int64) = c.value

function value(t::Table1D, iID::Int64)
    t.bound[] ||
        throw(ArgumentError("dependent value from $(t.source) is not bound to field \"$(t.field_name)\"; call bind_dependents! first"))
    x = t.field[][iID]
    if t.warn[] && (x < t.x[1] || x > t.x[end])
        @warn "$(t.field_name) = $x is outside the data range [$(t.x[1]), $(t.x[end])] of $(t.source). Using the nearest boundary value."
        t.warn[] = false
    end
    return evaluate(t.spline, x)
end

function bind_table!(t::Table1D, field::Vector{Float64})
    t.field[] = field
    t.bound[] = true
    return t
end

"""
    read_table(file, alias, path, ctx)

Reads the data file for parameter `alias`. Format: optional `#` comment lines,
a line `header: <field> <column> ...`, then whitespace separated numbers. The
first column is the node field the value depends on; the column used is
`alias` with spaces replaced by underscores. Problems are added to `ctx` at
`path` and `nothing` is returned.
"""
function read_table(file::String, alias::String, path::String, ctx::ParseContext)
    if !isfile(file)
        add_error!(ctx, path, "data file \"$file\" not found")
        return nothing
    end
    header = String[]
    rows = Vector{Vector{Float64}}()
    for (line_number, raw_line) in enumerate(eachline(file))
        # Windows editors may start the file with a UTF-8 byte order mark
        line = strip(line_number == 1 ? lstrip(==('﻿'), raw_line) : raw_line)
        (isempty(line) || startswith(line, "#")) && continue
        if isempty(header)
            if startswith(line, "header:")
                header = String.(split(line)[2:end])
                continue
            end
            add_error!(ctx, path,
                       "$file line $line_number: expected 'header: <field> <column> ...' before the data")
            return nothing
        end
        parts = split(line)
        if length(parts) != length(header)
            add_error!(ctx, path,
                       "$file line $line_number: expected $(length(header)) values, got $(length(parts))")
            return nothing
        end
        row = tryparse.(Float64, parts)
        if any(isnothing, row)
            add_error!(ctx, path, "$file line $line_number: non-numeric value")
            return nothing
        end
        push!(rows, Float64.(row))
    end
    column = replace(alias, " " => "_")
    index = findfirst(==(column), header)
    if index === nothing || index == 1
        add_error!(ctx, path,
                   "$file has no column \"$column\" (header: $(join(header, " ")))")
        return nothing
    end
    if length(rows) < 2
        add_error!(ctx, path, "$file needs at least 2 data rows")
        return nothing
    end
    x = [row[1] for row in rows]
    y = [row[index] for row in rows]
    if !all(diff(x) .> 0)
        add_error!(ctx, path,
                   "$file: first column ($(header[1])) must be strictly increasing")
        return nothing
    end
    return Table1D(header[1], x, y, file)
end

"""
    combine(f, a, b)

Applies `f` pointwise to two dependent values, e.g. to derive a shear modulus
from Young's modulus and Poisson's ratio in `derive`. Tables must depend on the
same field. The result is unbound.
"""
combine(f, a::Constant, b::Constant) = Constant(f(a.value, b.value))
combine(f, a::Table1D, b::Constant) = Table1D(a.field_name, a.x, f.(a.y, b.value), a.source)
combine(f, a::Constant, b::Table1D) = Table1D(b.field_name, b.x, f.(a.value, b.y), b.source)
function combine(f, a::Table1D, b::Table1D)
    a.field_name == b.field_name ||
        throw(ArgumentError("cannot combine values depending on \"$(a.field_name)\" ($(a.source)) and \"$(b.field_name)\" ($(b.source))"))
    x = sort!(unique!(vcat(a.x, b.x)))
    return Table1D(a.field_name, x, f.(evaluate(a.spline, x), evaluate(b.spline, x)),
                   "$(a.source) + $(b.source)")
end
