# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    ParamsDefinitionError(msg)

Raised when an `@params` declaration or a registration is invalid. This is a
programming error in a module, not a problem in the user's input deck.
"""
struct ParamsDefinitionError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ParamsDefinitionError) = print(io, "ParamsDefinitionError: ", e.msg)

"""
    InputError(path, message, severity)

One problem found in the input deck. `severity` is `:error` or `:warning`.
"""
struct InputError
    path::String
    message::String
    severity::Symbol
end

"""
    ParseContext(; directory = "", strict = true)

State shared while reading an input deck: the deck's directory (for relative
data file paths), strict mode, and all problems found so far.
"""
mutable struct ParseContext
    directory::String
    strict::Bool
    errors::Vector{InputError}
end
function ParseContext(; directory::AbstractString = "", strict::Bool = true)
    return ParseContext(String(directory), strict, InputError[])
end

function add_error!(ctx::ParseContext, path::AbstractString, msg::AbstractString)
    push!(ctx.errors, InputError(path, msg, :error))
end
function add_warning!(ctx::ParseContext, path::AbstractString, msg::AbstractString)
    push!(ctx.errors, InputError(path, msg, :warning))
end
has_errors(ctx::ParseContext) = any(e -> e.severity == :error, ctx.errors)

"""
    join_path(path, key)

Appends `key` to an error path. Keys that are not plain identifiers are quoted.
"""
function join_path(path::AbstractString, key::AbstractString)
    segment = occursin(r"^[A-Za-z0-9_]+$", key) ? String(key) : "\"$key\""
    return isempty(path) ? segment : "$path.$segment"
end

_normalize(s::AbstractString) = lowercase(filter(c -> isletter(c) || isdigit(c), s))

function levenshtein(a::AbstractString, b::AbstractString)
    a_chars, b_chars = collect(a), collect(b)
    n = length(b_chars)
    previous = collect(0:n)
    current = similar(previous)
    for (i, ca) in enumerate(a_chars)
        current[1] = i
        for j in 1:n
            cost = ca == b_chars[j] ? 0 : 1
            current[j + 1] = min(previous[j + 1] + 1, current[j] + 1, previous[j] + cost)
        end
        previous, current = current, previous
    end
    return previous[n + 1]
end

"""
    suggest(key, candidates)

Returns the candidate closest to `key` (ignoring case, spaces and
punctuation), or `nothing` if none is close enough.
"""
function suggest(key::AbstractString, candidates)
    normalized_key = _normalize(key)
    best = nothing
    best_distance = typemax(Int)
    for candidate in candidates
        distance = levenshtein(normalized_key, _normalize(candidate))
        if distance < best_distance
            best, best_distance = String(candidate), distance
        end
    end
    best === nothing && return nothing
    return best_distance <= max(2, length(normalized_key) ÷ 4) ? best : nothing
end

function unknown_key_message(key::AbstractString, candidates)
    suggestion = suggest(key, candidates)
    return suggestion === nothing ? "unknown key" :
           "unknown key — did you mean \"$suggestion\"?"
end

function format_errors(errors::Vector{InputError})
    io = IOBuffer()
    print(io, "Input errors (", length(errors), "):")
    for e in errors
        print(io, "\n  ", e.path, ": ", e.message)
    end
    if any(e -> startswith(e.message, "unknown key"), errors)
        print(io,
              "\nUnknown keys can be reported as warnings instead: run with --no_strict or set \"Strict Validation: false\".")
    end
    return String(take!(io))
end

"""
    report!(ctx)

Logs all warnings, then aborts with every error at once if there are any.
"""
function report!(ctx::ParseContext)
    for w in ctx.errors
        w.severity == :warning && @warn "$(w.path): $(w.message)"
    end
    errors = filter(e -> e.severity == :error, ctx.errors)
    isempty(errors) || @abort format_errors(errors)
    return nothing
end

"""
    strict_mode(input; no_strict_flag = false)

Strict validation is on unless the command line flag `--no_strict` is given or
the input deck sets `Strict Validation: false`.
"""
function strict_mode(input::AbstractDict; no_strict_flag::Bool = false)
    no_strict_flag && return false
    value = get(input, "Strict Validation", true)
    value isa Bool || @abort "\"Strict Validation\" must be true or false, got $(repr(value))"
    return value
end
