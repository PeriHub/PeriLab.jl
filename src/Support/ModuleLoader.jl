# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
# SPDX-FileCopyrightText: 2026
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Module discovery and loading: local, directory-scanned modules and
# license-server-fetched modules.
#
# INTEGRITY MODEL (for load_licensed_modules/licensed_modules) -- read this
# before changing it:
#   The download_token from /validate is signed with a server-only secret
#   (SIGNING_SECRET). This process has no way to verify that signature --
#   it can only present the token back to the server. That means the token
#   is NOT usable as a client-side integrity check on module content.
#   Instead, module content is verified against CONTENT_SIGNING_SECRET: a
#   SEPARATE secret this client deployment genuinely holds (set as a
#   runtime env var by whoever provisions this container, same as
#   LICENSE_SERVER_URL -- not the customer's own LICENSE_KEY). The
#   server signs each module's content with the same secret and returns the
#   signature in the X-Content-Signature response header; this loader
#   recomputes it locally and compares before ever calling include_string.
#
# USAGE -- two entry points for licensed modules, pick one per call site:
#   `licensed_modules(target_module, wanted_type)` is what feature modules
#   (Additive, a gcode/mesh-input module, ...) should call. `wanted_type`
#   is a category (e.g. "Additive", set per module in modules_config.yaml),
#   not a file name -- filtering by type means a single license feature can
#   bundle multiple files for different purposes (e.g. "additive" unlocks
#   both an Additive-typed module and a GcodeInput-typed one, routed to two
#   different parent modules) and a type can resolve to more than one file
#   (multiple modules sharing the same type all get loaded together).
#
#   Where the modules actually come from is controlled entirely by ENV, so
#   call sites never change between modes:
#     LOCAL DEV MODE -- if LICENSED_MODULES_CONFIG and LICENSED_MODULES_DIR
#     are both set, modules are read straight off local disk using the
#     SAME modules_config.yaml shape the server uses (features -> list of
#     {file, type} entries), with no license/feature gating at all -- if
#     you can see the file, you can load it, since this is for local
#     development. Point these at your own working copies of
#     server/modules_config.yaml and server/modules/ (or a checkout of
#     them) and edit a module's .jl file directly; changes are picked up
#     on your next Julia restart, no docker rebuild, no license server
#     round trip. Files are `include()`d (via `Base.include(target_module,
#     path)`), not `include_string`'d, so stack traces and debugging
#     tools point at the real file path. No content-signature check either
#     -- that's a network-transport concern and there's no network here.
#     If BOTH local and server env vars happen to be set, local mode wins.
#     Setting only one of LICENSED_MODULES_CONFIG/LICENSED_MODULES_DIR is
#     treated as a misconfiguration and throws.
#     SERVER MODE -- if LICENSE_SERVER_URL/LICENSE_KEY/
#     CONTENT_SIGNING_SECRET are set (and local mode isn't configured),
#     modules are fetched, verified, and include_string'd as described
#     under INTEGRITY MODEL above. The network round trip is cached for
#     the process lifetime (so N feature modules calling this costs one
#     validate+download, not N).
#     NEITHER configured: silent no-op, returns `Any[]`.
#
#   `load_licensed_modules(url, key, secret; target_module=...)` is the
#   low-level, explicit-argument, uncached, server-only version -- for
#   tests or advanced setups (e.g. more than one license server in the
#   same process). It has no local-mode equivalent.
#
# `find_registered_modules(directory, register_function)` finds the model
# files of a category by their `register_*` call; a missing directory has
# none.
#
# Both `licensed_modules` and `find_module_files` return a `Vector{Any}` of
# `Dict("File" => ..., "Module Name" => ...)` entries in the shape
# `create_module_specifics` expects, so their results concatenate directly.

module ModuleLoader

using ..PeriLabExceptions: @abort
using DotEnv
DotEnv.load!()

import Downloads
import YAML
import JSON3
import SHA

export find_module_files, find_registered_modules, create_module_specifics,
       licensed_modules

"""
	find_jl_files(directory::AbstractString)

Recursively find Julia files (.jl) in a directory.

This function recursively searches for Julia source files with the ".jl" extension
in the specified directory and its subdirectories. It returns a vector of file paths
for all the found .jl files.

# Arguments
- `directory::AbstractString`: The directory in which to search for .jl files.

# Returns
A vector of strings, where each string is a file path to a .jl file found in the
specified directory and its subdirectories.

# Example
```julia
jl_files = find_jl_files("/path/to/modules")
for jl_file in jl_files
	println("Found Julia file: ", jl_file)
end
```
"""
function find_jl_files(directory::AbstractString)
    if !isdir(directory)
        @abort "$directory does not exists. Modules won't be loaded accurately."
    end
    return sort!([joinpath(root, file)
                  for (root, _, files) in walkdir(directory)
                  for file in files
                  if endswith(file, ".jl")])
end

"""
	find_module_files(directory::AbstractString, specific::String)

Search for Julia modules containing a specific function in a given directory.

This function searches for Julia modules (files with `.jl` extension) in the specified
directory and checks if they contain a specific function. It returns a list of dictionaries
where each dictionary contains the file path and the name of the module where the specific
function is found.

# Arguments
- `directory::AbstractString`: The directory to search for Julia modules.
- `specific::String`: The name of the specific function to search for.

# Returns
An array of dictionaries, where each dictionary has the following keys:
- `"File"`: The file path to the module where the specific function is found.
- `"Module Name"`: The name of the module where the specific function is found.

# Example
```julia
result = find_module_files("/path/to/modules", "my_function")
for module_info in result
	println("Function found in module: ", module_info["Module Name"])
	println("Module file path: ", module_info["File"])
end
```
"""
function find_module_files(directory::AbstractString, specific::String)
    files_in_folder = find_jl_files(directory)
    module_list = []
    module_name = ""
    for filename in files_in_folder
        file = open(filename, "r")
        for line in eachline(file)
            if occursin(r"\bmodule\b", line)
                module_name = split(line)[2]
            end
            if occursin("function " * specific * "()", line)
                push!(module_list, Dict("File" => filename, "Module Name" => module_name))
                break
            end
        end
        close(file)
    end
    return module_list
end

"""
    find_registered_modules(directory, register_function; exclude = String[])

The model files below `directory` that call `register_function` (e.g.
`"register_damage"`), the call in a comment included, so that the unregistered
templates load as well. Files and directories in `exclude` are skipped; a missing
`directory` has no model files.
"""
function find_registered_modules(directory::AbstractString, register_function::String;
                                 exclude = String[])
    isdir(directory) || return String[]
    call = register_function * "("
    excluded(file) = any(e -> file == e || startswith(file, joinpath(e, "")), exclude)
    return [file
            for file in find_jl_files(directory)
            if !excluded(file) && any(line -> occursin(call, line), eachline(file))]
end

"""
	create_module_specifics(name::String, module_list::Dict{String,AbstractString}(),specifics::Dict{String,String}(), values::Tuple)

Searches for a specific function within a list of modules and calls that function if found.

This function iterates over a list of modules specified in `module_list` and looks for a module-specific function specified in the `specifics` dictionary. If the module and function are found, it calls that function with the provided `values` tuple.

# Arguments
- `name::String`: The name to match against the module names.
- `module_list::Dict{String, AbstractString}`: A dictionary of module names mapped to abstract strings.
- `specifics::Dict{String, String}`: A dictionary specifying the module-specific function to call for each module.
- `values::Tuple`: A tuple of values to be passed as arguments to the module-specific function.

# Example
```julia
module_list = Dict("Module1" => "Module1Name", "Module2" => "Module2Name")
specifics = Dict("Module1Name" => "module1_function", "Module2Name" => "module2_function")
values = (arg1, arg2)
create_module_specifics("Module1Name", module_list, specifics, values)
```
"""
function create_module_specifics(name::Union{String,SubString},
                                 module_list::Vector{Any},
                                 own_module::Module,
                                 specifics::Dict{String,String},
                                 values::Tuple)
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * specifics["Name"] * "()"
        if Base.eval(own_module, Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"] * "." * specifics["Call Function"]
            function_call = Base.eval(own_module, Meta.parse(parse_statement))
            return function_call(values...)
        end
    end
    return nothing
end
"""
	create_module_specifics(name::String, module_list::Dict{String,AbstractString}(),specifics::Dict{String,String}())
	# Returns: the function itself
"""
function create_module_specifics(name::Union{String,SubString},
                                 module_list::Vector{Any},
                                 own_module::Module,
                                 specifics::Dict{String,String})
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * specifics["Name"] * "()"
        if Base.eval(own_module, Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"] * "." * specifics["Call Function"]
            function_call = Base.eval(own_module, Meta.parse(parse_statement))
            return function_call
        end
    end
    return nothing
end
# only module
function create_module_specifics(name::Union{String,SubString},
                                 module_list::Vector{Any},
                                 own_module::Module,
                                 get_model_name::String)
    for m in module_list
        parse_statement = "module_name=" * m["Module Name"] * "." * get_model_name * "()"
        if Base.eval(own_module, Meta.parse(parse_statement)) == name
            parse_statement = m["Module Name"]
            module_call = Base.eval(own_module, Meta.parse(parse_statement))
            return module_call
        end
    end
    return nothing
end

# Process-wide cache for license-server-fetched modules. The network round trip (validate + download + verify) doesn't depend on
# which feature module is asking, so it happens at most once per process.
# What DOES depend on the caller is which Julia namespace the source gets
# include_string'd into -- create_module_specifics later does
# Base.eval(own_module, ...), so the module must exist inside whichever
# module actually calls licensed_modules. _INCLUDED tracks that per
# (target_module, module_name) pair so a repeat call with the same target
# doesn't re-include (and spam a "replacing module" warning).
const _Entry = NamedTuple{(:module_name, :module_type, :source, :detected_name),
                          Tuple{String,String,String,String}}
const _FETCH_CACHE = Ref{Union{Nothing,Vector{_Entry}}}(nothing)
const _INCLUDED = IdDict{Module,Set{String}}()
const _CACHE_LOCK = ReentrantLock()
const _WARNED_NOT_CONFIGURED = Ref(false)

# """
#     __init__()

# Resets all cache state fresh on every new process. This matters
# specifically because precompilation (a plain `Pkg.precompile()`, an app
# registration step like `Pkg.Apps.develop`, or a PackageCompiler sysimage
# build) can execute a package's real entry point once as part of building
# it -- if `licensed_modules()` happens to run during that step (e.g.
# because it runs before any license env vars exist, during `docker
# build`), the resulting "no license configured" state would otherwise get
# frozen into the precompiled cache and never be re-checked again at actual
# runtime, no matter what env vars a container is started or exec'd with
# later.

# `__init__()` runs fresh every time a module is loaded into a NEW Julia
# process, regardless of what happened during precompilation -- that's
# precisely what it's for, and it doesn't require discarding or
# regenerating any compiled code to get a correct, fresh check on each real
# invocation. This is the reason a `.ji`-cache-deletion workaround in a
# container entrypoint shouldn't be necessary: it forces a full recompile to
# get a fresh check, when only the cache in this module needs resetting.

"""
    licensed_modules(target_module::Module, wanted_type=nothing;
                      force_refresh::Bool=false, machine_id=gethostname())

Reads `LICENSE_SERVER_URL`, `LICENSE_KEY`, and
`CONTENT_SIGNING_SECRET` from `ENV` and `include_string`s licensed modules
into `target_module`. This is the entry point feature modules (Additive,
a gcode/mesh-input module, ...) should call -- no need to redefine the
URL/key/secret in every file.

`wanted_type`, if given, is a category tag (e.g. `"Additive"`), not a file
name -- the server assigns a `type` to every module in
`modules_config.yaml`, and `/validate` returns that mapping. **This
matters for two reasons**: first, a single license feature can unlock more
than one file for different purposes -- e.g. the `additive` feature
unlocks both `AdditiveMaterialModel.jl` (type `"Additive"`) and
`GcodeMeshInput.jl` (type `"GcodeInput"`) together, but those two files
belong in two different PeriLab parent modules. Second, more than one file
can share the same type -- e.g. two different additive material model
implementations both tagged `"Additive"` -- and a caller asking for that
type gets all of them, not just one. Without `wanted_type`, EVERY module
the license entitles gets `include_string`'d into EVERY call site that
asks for `target_module`. Passing `wanted_type` filters that: only modules
tagged with it are `include_string`'d into `target_module`, and only they
appear in the returned list. A type the license doesn't grant anything for
resolves to zero modules (not an error) -- same "optional" philosophy as
the rest of this module -- with a `@debug` note, to help you notice a
typo'd type versus a genuinely unlicensed one.

The underlying validate+download+verify round trip is unaffected by
`wanted_type` and still runs at most once per process regardless of how
many call sites there are, each asking for their own type -- filtering
happens only at the `include_string` step, not the network fetch.

- `LICENSE_KEY` not set: silently returns `Any[]`. No license key means no
  licensed modules, not an error -- this is what makes licensing fully
  optional. `LICENSE_SERVER_URL` defaults to the PeriHub license server.
- `LICENSE_KEY` set but `CONTENT_SIGNING_SECRET` missing: that's a
  misconfiguration, not "no license" -- throws rather than silently
  disabling the integrity check.
- Both set: the actual validate+download+verify round trip runs at
  most once per process no matter how many call sites there are.

Set `force_refresh=true` to bypass the cache (e.g. to pick up a renewed
license without restarting the process).
"""
function licensed_modules(target_module::Module,
                          wanted_type::Union{Nothing,AbstractString} = nothing;
                          force_refresh::Bool = false,
                          machine_id::AbstractString = gethostname())
    DotEnv.load!()

    local_config = get(ENV, "LICENSED_MODULES_CONFIG", "")
    local_dir = get(ENV, "LICENSED_MODULES_DIR", "")

    if !isempty(local_config) || !isempty(local_dir)
        if isempty(local_config) || isempty(local_dir)
            @abort "Only one of LICENSED_MODULES_CONFIG/LICENSED_MODULES_DIR is set -- " *
                   "set both to enable local dev mode, or neither to use the license server."
        end
        return _licensed_modules_local(target_module, wanted_type, local_config, local_dir)
    end

    url = get(ENV, "LICENSE_SERVER_URL", "https://perihub-license.nimbus-extern.dlr.de")
    key = get(ENV, "LICENSE_KEY", "")
    secret = get(ENV, "CONTENT_SIGNING_SECRET", "")

    if isempty(key)
        lock(_CACHE_LOCK) do
            if !_WARNED_NOT_CONFIGURED[]
                @info "LICENSE_KEY not set -- skipping licensed modules"
                _WARNED_NOT_CONFIGURED[] = true
            end
        end
        return Any[]
    end

    if isempty(secret)
        @abort "LICENSE_KEY is set but " *
               "CONTENT_SIGNING_SECRET is not -- refusing to load licensed " *
               "modules without a way to verify their integrity."
    end

    entries = lock(_CACHE_LOCK) do
        if force_refresh || _FETCH_CACHE[] === nothing
            _FETCH_CACHE[] = _fetch_and_verify_all(url, key, secret, machine_id)
        end
        return _FETCH_CACHE[]
    end

    module_list = Vector{Any}()
    matched_any = false
    lock(_CACHE_LOCK) do
        included = get!(_INCLUDED, target_module, Set{String}())
        for entry in entries
            if wanted_type !== nothing && entry.module_type != wanted_type
                continue
            end
            matched_any = true
            if entry.module_name ∉ included
                # Defines the module directly in target_module's namespace.
                # Using the remote module_name as the "filename" argument
                # only affects stack traces/error messages -- it is not
                # read from or written to disk.
                Base.include_string(target_module, entry.source, entry.module_name)
                push!(included, entry.module_name)
            end
            push!(module_list,
                  Dict("File" => "license://" * entry.module_name,
                       "Module Name" => entry.detected_name))
        end
    end

    if wanted_type !== nothing && !matched_any
        @debug "No licensed modules of the requested type are available to this license" wanted_type
    end

    return module_list
end

"""
    _licensed_modules_local(target_module, wanted_type, config_path, modules_dir)

Local dev mode backing `licensed_modules` when `LICENSED_MODULES_CONFIG`
and `LICENSED_MODULES_DIR` are both set. Reads the same
`features -> modules: [{file, type}]` shape the license server's
`modules_config.yaml` uses, but applies NO license/feature gating -- every
module in the config that matches `wanted_type` (or every module, if
`wanted_type` is `nothing`) is loaded. This is intentional: local mode is
for developing/testing licensed modules on your own machine, where you
already have full access to the files.

Re-reads the config and files fresh on every call (no caching, unlike
server mode) -- cheap since it's local disk, and means edits are visible
without restarting anything except the Julia session that already
`include()`d the previous version (Julia doesn't support redefining an
already-loaded module in place any more than server mode does).

Files are loaded with `Base.include(target_module, path)`, not
`include_string`, specifically so error messages and debugger/backtrace
tooling show the real file path during local development.
"""
function _licensed_modules_local(target_module::Module,
                                 wanted_type::Union{Nothing,AbstractString},
                                 config_path::AbstractString,
                                 modules_dir::AbstractString)
    if !isfile(config_path)
        @abort "LICENSED_MODULES_CONFIG is set to '$(config_path)' but that file doesn't exist"
    end
    if !isdir(modules_dir)
        @abort "LICENSED_MODULES_DIR is set to '$(modules_dir)' but that directory doesn't exist"
    end

    data = YAML.load_file(config_path)
    features = get(data, "features", Dict())

    included = lock(_CACHE_LOCK) do
        get!(_INCLUDED, target_module, Set{String}())
    end

    module_list = Vector{Any}()
    seen = Set{String}()
    matched_any = false

    for (_, feature) in features
        for m in get(feature, "modules", Any[])
            file_name = m["file"]
            file_name in seen && continue
            push!(seen, file_name)

            module_type = string(get(m, "type", file_name))
            if wanted_type !== nothing && module_type != wanted_type
                continue
            end
            matched_any = true

            path = joinpath(modules_dir, file_name)
            if !isfile(path)
                @warn "Local licensed module config references '$(file_name)' but it doesn't exist at $(path) -- skipping" path
                continue
            end

            if file_name ∉ included
                Base.include(target_module, path)
                push!(included, file_name)
            end

            detected_name = _detect_module_name(read(path, String))
            push!(module_list, Dict("File" => path, "Module Name" => detected_name))
        end
    end

    if wanted_type !== nothing && !matched_any
        @debug "No local licensed modules of the requested type are configured" wanted_type config_path
    end

    return module_list
end

function _fetch_and_verify_all(license_server_url::AbstractString,
                               license_key::AbstractString,
                               content_signing_secret::AbstractString,
                               machine_id::AbstractString)
    validation = _validate(license_server_url, license_key, machine_id)

    if !validation.valid
        @abort "PeriHub license check failed: $(get(validation, :reason, "unknown reason"))"
    end

    # Module downloads use the short-lived, module-scoped download_token
    # returned by /validate -- not the long-lived license_key itself -- so
    # the raw key is sent over the wire only once per validation.
    download_token = validation.download_token

    entries = _Entry[]
    for module_name in validation.modules
        source,
        signature = _fetch_module_source(license_server_url, download_token, module_name)

        expected_signature = bytes2hex(SHA.hmac_sha256(Vector{UInt8}(content_signing_secret),
                                                       source))
        if isempty(signature)
            @abort "Server did not return X-Content-Signature for module '$(module_name)' -- refusing to load it"
        elseif signature != expected_signature
            @abort "Content signature mismatch for module '$(module_name)': " *
                   "expected $(expected_signature), got $(signature). Refusing to " *
                   "load it (module may have been tampered with in transit)."
        end

        # Falls back to the file name itself as its own type if the server
        # config didn't tag one -- defensive, not expected in normal use.
        module_type = string(get(validation.module_types, Symbol(module_name), module_name))
        detected_name = _detect_module_name(source)
        push!(entries,
              (module_name = module_name, module_type = module_type,
               source = source, detected_name = detected_name))
    end
    return entries
end

# HTTP request; errors on a non-2xx status. Returns the response and its body.
function _request(url::AbstractString; method = "GET", headers = Pair{String,String}[],
                  body = nothing)
    output = IOBuffer()
    resp = Downloads.request(url; method = method, headers = headers, output = output,
                             input = body === nothing ? nothing : IOBuffer(body))
    200 <= resp.status < 300 ||
        @abort "$method $url failed with HTTP status $(resp.status)"
    return resp, take!(output)
end

function _validate(license_server_url::AbstractString, license_key::AbstractString,
                   machine_id::AbstractString)
    body = JSON3.write(Dict("license_key" => license_key, "machine_id" => machine_id))
    _, response = _request(rstrip(license_server_url, '/') * "/api/v1/validate";
                           method = "POST",
                           headers = ["Content-Type" => "application/json"], body = body)
    return JSON3.read(response)
end

function _fetch_module_source(license_server_url::AbstractString,
                              download_token::AbstractString,
                              module_name::AbstractString)
    resp, source = _request(rstrip(license_server_url, '/') * "/api/v1/modules/" *
                            module_name;
                            headers = ["Authorization" => "Bearer " * download_token])
    i = findfirst(h -> lowercase(first(h)) == "x-content-signature", resp.headers)
    signature = i === nothing ? "" : last(resp.headers[i])
    return String(source), signature
end

function _detect_module_name(source::AbstractString)
    for line in split(source, '\n')
        if occursin(r"\bmodule\b", line)
            return split(line)[2]
        end
    end
    @abort "Downloaded module source has no top-level `module ... end` block"
end

end # module ModuleLoader
