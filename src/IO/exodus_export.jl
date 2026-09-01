# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus
using Dates
using DataStructures

export get_paraview_coordinates
export create_result_file
export init_results_in_exodus
export merge_exodus_file
export compute_bond_connectivity
export bond_export_sizes
export init_bond_information_export
export resolve_block_selection
export check_block_selection
export element_block_sizes
export write_bond_results_in_exodus
export BondBlock

"""
    BondBlock

Holds everything needed to write one bond block, built once in
[`compute_bond_connectivity`](@ref).

The three fields are aligned column by column: bond `k` runs from `nodes[k]` to the
second entry of `conn[:, k]`, and corresponds to entry `neighbor_index[k]` of the
owning node's neighbourhood list. Keeping the neighbour index around is what makes the
result output reproduce exactly the order used for the connectivity — otherwise bond
values and bond elements silently drift apart.

Bonds are directed, so the same node pair may appear twice, once per direction, each
with its own values. The owning node is always `conn[1, k]`.

# Fields
- `conn::Matrix{Int64}`: `2 x n_bonds` BAR2 connectivity
- `nodes::Vector{Int64}`: owning node of each bond
- `neighbor_index::Vector{Int64}`: position of the bond inside `nlist[node]`
"""
struct BondBlock
    conn::Matrix{Int64}
    nodes::Vector{Int64}
    neighbor_index::Vector{Int64}
end

Base.length(b::BondBlock) = length(b.nodes)

"""
    create_result_file(filename, num_nodes, num_dim, num_elem_blks, num_node_sets,
                       num_elements = 0, FEtopology = nothing)

Creates an exodus file for the results.

!!! important
    `num_elem_blks` must be the **total** number of element blocks that will later be
    written, i.e. the number of material blocks **plus** the number of non-empty bond
    blocks. Use [`bond_export_sizes`](@ref) to obtain the latter. exodusII fixes this
    count at initialisation; writing more blocks than declared makes `ex_put_block`
    fail.

# Arguments
- `filename::Union{AbstractString,String}`: The name of the file to create
- `num_nodes::Int64`: The number of nodes
- `num_dim::Int64`: The number of dimensions
- `num_elem_blks::Int64`: Total number of element blocks (material + bond blocks)
- `num_node_sets::Int64`: The number of node sets
- `num_elements::Int64`: Number of additional elements (FE elements + bond elements)
- `FEtopology::Union{Nothing,Matrix{Int64}}`: FE topology, if a FE part is present
# Returns
- `result_file::Dict{String,Any}`: A dictionary containing the filename and the exodus file
"""
function create_result_file(filename::Union{AbstractString,String},
                            num_nodes::Int64,
                            num_dim::Int64,
                            num_elem_blks::Int64,
                            num_node_sets::Int64,
                            num_elements::Int64 = 0,
                            FEtopology::Union{Nothing,Matrix{Int64}} = nothing)
    if isfile(filename)
        rm(filename)
    end
    maps_int_type = Int32
    ids_int_type = Int32
    bulk_int_type = Int32
    float_type = Float64

    num_elems = num_nodes + num_elements
    if !isnothing(FEtopology)
        element_nodes = unique(reduce(vcat, FEtopology))
        num_elems -= length(element_nodes)
    end

    num_side_sets = 0

    init = Initialization{Int32(num_dim),
                          Int32(num_nodes),
                          Int32(num_elems),
                          Int32(num_elem_blks),
                          Int32(num_node_sets),
                          Int32(num_side_sets)}()

    @info "Create output " * filename
    @debug "Exodus init: $num_nodes nodes, $num_elems elements, $num_elem_blks blocks"

    exo_db = ExodusDatabase{maps_int_type,ids_int_type,bulk_int_type,float_type}(filename,
                                                                                 "w",
                                                                                 init)
    return Dict("filename" => filename, "file" => exo_db, "type" => "Exodus")
end

"""
    paraview_specifics(dof::Int64)

Returns the paraview specific dof

# Arguments
- `dof::Int64`: The degrees of freedom
# Returns
- `paraview_specifics::String`: The paraview specific dof
"""
function paraview_specifics(dof::Int64)
    convention = Dict(1 => "x", 2 => "y", 3 => "z")
    return convention[dof]
end

"""
    get_paraview_coordinates(dof::Int64, refDof::Int64)

Returns the paraview specific dof

# Arguments
- `dof::Int64`: The degrees of freedom
- `refDof::Int64`: The reference degrees of freedom
# Returns
- `paraview_specifics::String`: The paraview specific dof
"""
function get_paraview_coordinates(dof::Int64, refDof::Int64)
    if dof > refDof
        @warn "Reference dof are to small and set to used dof"
        refDof = dof
    end
    if refDof < 4
        return paraview_specifics(dof)
    end
    if refDof < 10
        return paraview_specifics(Int(ceil(dof / 3))) *
               paraview_specifics(dof - Int(ceil(dof / 3 - 1) * 3))
    end

    @abort "not yet exportable as one variable"
end

"""
    get_block_nodes(block_Id::AbstractVector{Int64}, block::Int64)

Returns the nodes of a block as a `1 x n` matrix, which is the layout exodusII
expects for a SPHERE block (one node per element).

# Arguments
- `block_Id::AbstractVector{Int64}`: The block Id
- `block::Int64`: The block
# Returns
- `nodes::Matrix{Int64}`: The nodes of the block, shaped `1 x n`
"""
function get_block_nodes(block_Id::AbstractVector{Int64}, block::Int64)
    conn = findall(x -> x == block, block_Id)
    return reshape(conn, 1, length(conn))
end

"""
    compute_bond_connectivity(block_Id, n_blocks, nlist; blocks = nothing)

Builds the BAR2 connectivity for the bond export **once**, so that the element count
used at initialisation, the connectivity written into the file, and the result values
written per time step all follow the same ordering.

Bonds are **directed**, which is what the peridynamic data structure actually holds:
`nlist[i]` containing `j` does not imply that `nlist[j]` contains `i`, and even when
both exist they carry separate state (bond damage, bond forces and so on differ per
direction). Each entry of each neighbourhood list therefore becomes its own BAR2
element, owned by the block of node `i`. Two elements on the same node pair are legal
in exodusII; they simply overlap geometrically.

Neighbours pointing outside the local node range (ghost entries of a decomposed mesh)
are skipped, because exodusII validates connectivity against the node count of the
file and those nodes carry no coordinates here.

# Arguments
- `block_Id::AbstractVector{Int64}`: Block id per node
- `n_blocks::Int64`: Number of material blocks
- `nlist::AbstractVector`: Neighbourhood list
- `blocks`: Optional subset of block indices to export bonds for
# Returns
- `bond_blocks::OrderedDict{Int64,BondBlock}`: Block index -> bond block.
  Blocks without bonds are absent from the dictionary.
"""
function compute_bond_connectivity(block_Id::AbstractVector{Int64},
                                   n_blocks::Int64,
                                   nlist::AbstractVector;
                                   blocks::Union{Nothing,AbstractVector{Int64}} = nothing)
    bond_blocks = OrderedDict{Int64,BondBlock}()
    n_nodes = length(block_Id)
    n_skipped = 0

    for block in 1:n_blocks
        if !isnothing(blocks) && !(block in blocks)
            continue
        end

        conn_nodes = get_block_nodes(block_Id, block)

        owners = Int64[]
        neighbor_indices = Int64[]
        partners = Int64[]

        for i in conn_nodes
            i > length(nlist) && continue
            neighbors = nlist[i]
            for (m, j) in enumerate(neighbors)
                if j > n_nodes
                    # ghost node of a decomposed mesh, not present in this file
                    n_skipped += 1
                    continue
                end
                push!(owners, i)
                push!(neighbor_indices, m)
                push!(partners, j)
            end
        end

        isempty(owners) && continue

        conn = Matrix{Int64}(undef, 2, length(owners))
        conn[1, :] .= owners
        conn[2, :] .= partners

        bond_blocks[block] = BondBlock(conn, owners, neighbor_indices)
    end

    if n_skipped > 0
        @debug "Bond export: $n_skipped bond(s) point to nodes outside this rank and " *
               "are not written."
    end

    return bond_blocks
end

"""
    bond_export_sizes(bond_blocks)

Returns the numbers needed by [`create_result_file`](@ref) for a given bond
connectivity.

# Arguments
- `bond_blocks`: Result of [`compute_bond_connectivity`](@ref)
# Returns
- `n_bond_elements::Int64`: Total number of BAR2 elements
- `n_bond_blocks::Int64`: Number of additional element blocks
"""
function bond_export_sizes(bond_blocks::AbstractDict{Int64,BondBlock})
    n_bond_elements = sum(length(b) for b in values(bond_blocks); init = 0)
    return n_bond_elements, length(bond_blocks)
end

"""
    resolve_block_selection(blocks)

Normalises the `"Bond Blocks"` entry into a list of block indices.

The entry is written as a plain, whitespace separated list of indices in the input
deck, e.g. `Bond Blocks: 1 3`. A single index and an already parsed collection of
integers are accepted as well. `nothing` or an empty entry means "all blocks".

# Arguments
- `blocks`: The selection taken from the output parameters
# Returns
- `Union{Nothing,Vector{Int64}}`: Block indices, or `nothing` for "all blocks"
"""
function resolve_block_selection(blocks)
    isnothing(blocks) && return nothing

    if blocks isa Integer
        return Int64[blocks]
    end

    if blocks isa AbstractString
        # "1 3" and "1, 3" both work; commas are treated as separators.
        tokens = split(replace(String(blocks), ',' => ' '))
        isempty(tokens) && return nothing
        indices = Int64[]
        for token in tokens
            index = tryparse(Int64, token)
            if isnothing(index)
                @abort "\"Bond Blocks\" expects block indices separated by spaces, " *
                       "got '$token'."
            end
            push!(indices, index)
        end
        return indices
    end

    isempty(blocks) && return nothing
    return collect(Int64, blocks)
end

"""
    check_block_selection(blocks, n_blocks)

Verifies that every index of a `"Bond Blocks"` selection refers to an existing block.

[`resolve_block_selection`](@ref) only parses; it has no idea which blocks exist.
Without this check an out-of-range index is silently skipped by the export loop, and
the user gets fewer bond blocks than asked for with no hint why.

# Arguments
- `blocks::Union{Nothing,AbstractVector{Int64}}`: The parsed selection, `nothing` for all blocks
- `n_blocks::Int64`: Number of existing material blocks
# Returns
- `blocks`: The unchanged input, so the call can be chained
"""
function check_block_selection(blocks, n_blocks::Int64)
    isnothing(blocks) && return blocks

    unknown = filter(b -> b < 1 || b > n_blocks, blocks)
    if !isempty(unknown)
        @abort "\"Bond Blocks\" references block(s) $(join(unknown, ", ")), " *
               "but only blocks 1 to $n_blocks exist."
    end

    return blocks
end

"""
    init_bond_information_export(block_Id, n_blocks, nlist, parameter)

Convenience wrapper: evaluates the `"Bond Export"` parameter and returns the bond
blocks together with the sizes required at initialisation.

# Arguments
- `block_Id::AbstractVector{Int64}`: Block id per node
- `n_blocks::Int64`: Number of material blocks
- `nlist::AbstractVector`: Neighbourhood list
- `parameter::Dict`: Output parameter block, may contain `"Bond Export"` and
  `"Bond Blocks"` (whitespace separated indices, e.g. `1 3`). Without `"Bond Blocks"`
  every block is exported.
# Returns
- `bond_blocks::OrderedDict{Int64,BondBlock}`: Bond blocks (empty if disabled)
- `n_bond_elements::Int64`: Total number of BAR2 elements
- `n_bond_blocks::Int64`: Number of additional element blocks
"""
function init_bond_information_export(block_Id::AbstractVector{Int64},
                                      n_blocks::Int64,
                                      nlist::AbstractVector,
                                      parameter::Dict)
    if !get(parameter, "Bond Export", false)
        return OrderedDict{Int64,BondBlock}(), 0, 0
    end

    blocks = resolve_block_selection(get(parameter, "Bond Blocks", nothing))
    check_block_selection(blocks, n_blocks)
    bond_blocks = compute_bond_connectivity(block_Id, n_blocks, nlist; blocks = blocks)
    n_bond_elements, n_bond_blocks = bond_export_sizes(bond_blocks)

    @info "Exporting bond information to exodus file: " *
          "$n_bond_elements bonds in $n_bond_blocks block(s)"

    return bond_blocks, n_bond_elements, n_bond_blocks
end

"""
    element_block_sizes(block_Id, n_blocks, bond_blocks)

Returns the number of elements per exodus block id, for both the material (SPHERE)
blocks and the bond (BAR2) blocks.

Exodus.jl exposes no truth table API, so an element variable has to be written for
**every** element block, not only for the bond blocks. This mapping supplies the
lengths of the zero vectors used for the material blocks.

# Arguments
- `block_Id::AbstractVector{Int64}`: Block id per node
- `n_blocks::Int64`: Number of material blocks
- `bond_blocks::AbstractDict{Int64,BondBlock}`: The bond blocks
# Returns
- `sizes::OrderedDict{Int64,Int64}`: exodus block id -> number of elements
"""
function element_block_sizes(block_Id::AbstractVector{Int64},
                             n_blocks::Int64,
                             bond_blocks::AbstractDict{Int64,BondBlock})
    sizes = OrderedDict{Int64,Int64}()
    for block in 1:n_blocks
        n = count(==(block), block_Id)
        n == 0 && continue
        sizes[block] = n
    end
    for (block, bond_block) in bond_blocks
        sizes[n_blocks + block] = length(bond_block)
    end
    return sizes
end

"""
    gather_bond_values(field, bond_block, dof)

Extracts the values of one bond variable in exactly the order of the bond
connectivity.

Two storage layouts are supported, decided once per field from a sample entry:

- **nested** (`Vector{Vector{...}}`, what PeriLab uses for bond fields):
  `field[node][neighbor]` is the value itself for a scalar field, a `Vector` for a
  vector field such as `Deformed Bond Geometry`, and a `Matrix` for a tensor field.
- **flat** (a `Matrix` or higher array per node): `field[node][neighbor, ...]`.

# Arguments
- `field`: The bond field
- `bond_block::BondBlock`: The bond block
- `dof`: `nothing` for scalars, an `Integer` for vectors, a 2-element collection for tensors
# Returns
- `values::Vector{Float64}`: One value per bond
"""
function gather_bond_values(field, bond_block::BondBlock, dof)
    n = length(bond_block)
    values = zeros(Float64, n)
    n == 0 && return values

    # Decide the layout once instead of branching inside the loop.
    nested = ndims(field[bond_block.nodes[1]]) == 1

    if nested
        @inbounds for k in 1:n
            entry = field[bond_block.nodes[k]][bond_block.neighbor_index[k]]
            if isnothing(dof)
                values[k] = Float64(entry)
            elseif dof isa Integer
                values[k] = Float64(entry[dof])
            else
                values[k] = Float64(entry[dof[1], dof[2]])
            end
        end
    else
        @inbounds for k in 1:n
            node_entry = field[bond_block.nodes[k]]
            m = bond_block.neighbor_index[k]
            if isnothing(dof)
                values[k] = Float64(node_entry[m])
            elseif dof isa Integer
                values[k] = Float64(node_entry[m, dof])
            else
                values[k] = Float64(node_entry[m, dof[1], dof[2]])
            end
        end
    end

    return values
end

"""
    init_results_in_exodus(exo, dof, output, coords, block_Id, all_block_name_list,
                           nsets, global_ids, PERILAB_VERSION, qa_vector;
                           bond_blocks = ..., bond_output_names = String[],
                           fem_block = nothing, topology = nothing,
                           elem_global_ids = nothing)

Initializes the results in exodus.

!!! note "Changed interface"
    The former positional arguments `bond_export::Bool` and `nlist` are gone. Pass the
    pre-computed `bond_blocks` from [`init_bond_information_export`](@ref) instead — the
    same object whose sizes were used for `create_result_file`. The remaining optional
    arguments are keywords now, which also removes the illegal ordering of the old
    signature (an optional argument may not precede a required one).

# Arguments
- `exo::ExodusDatabase`: The exodus database
- `dof`: Degrees of freedom
- `output::Dict`: The output definition
- `coords::Union{Matrix{Int64},Matrix{Float64}}`: The coordinates
- `block_Id::Vector{Int64}`: The block Id per node
- `all_block_name_list::Vector{String}`: The block names
- `nsets::Dict{String,Vector{Int64}}`: The node sets
- `global_ids::Vector{Int64}`: The global ids
- `PERILAB_VERSION::String`: PeriLab version string
- `qa_vector::Vector{String}`: Additional QA records (at most 2 are stored)
# Keywords
- `bond_blocks`: Bond blocks, see [`compute_bond_connectivity`](@ref)
- `bond_output_names::Vector{String}`: Names of the bond variables, written as element variables
- `fem_block::Union{Nothing,Vector{Bool}}`: Per-node flag marking FE nodes
- `topology::Union{Nothing,Matrix{Int64}}`: FE topology
- `elem_global_ids::Union{Nothing,Vector{Int64}}`: Global element ids for the FE part
# Returns
- `exo::ExodusDatabase`: The exodus file
"""
function init_results_in_exodus(exo::ExodusDatabase,
                                dof,
                                output::Dict,
                                coords::Union{Matrix{Int64},Matrix{Float64}},
                                block_Id::Vector{Int64},
                                all_block_name_list::Vector{String},
                                nsets::Dict{String,Vector{Int64}},
                                global_ids::Vector{Int64},
                                PERILAB_VERSION::String,
                                qa_vector::Vector{String};
                                bond_blocks::AbstractDict{Int64,BondBlock} = OrderedDict{Int64,
                                                                                         BondBlock}(),
                                bond_output_names::Vector{String} = String[],
                                fem_block::Union{Nothing,Vector{Bool}} = nothing,
                                topology::Union{Nothing,Matrix{Int64}} = nothing,
                                elem_global_ids::Union{Nothing,Vector{Int64}} = nothing)
    qa = Matrix{String}(undef, 1, 4)
    # Only 4 entries with 32 chars possible!
    qa[1] = "PeriLab $PERILAB_VERSION"
    qa[2] = Dates.format(Dates.now(), "mm/dd/yyyy HH:MM")
    if length(qa_vector) > 2
        @warn "Only the first two QA entries are stored, $(length(qa_vector)) were given."
    end
    for (id, entry) in enumerate(qa_vector)
        id > 2 && break
        qa[2 + id] = entry
    end
    write_qa(exo, qa)

    info = [
        "PeriLab Version $PERILAB_VERSION, under BSD License",
        "Copyright (c) 2023, Christian Willberg, Jan-Timo Hesse",
        "compiled with Julia Version " * string(VERSION)
    ]
    write_info(exo, info)

    # check if type of coords is int or Float64
    if typeof(coords) in [Matrix{Int64}, Matrix{Float64}]
        coords = convert(Matrix{Float64}, coords)
    end

    write_coordinates(exo, coords)

    id::Int32 = 0
    for name in eachindex(nsets)
        id += Int32(1)
        nsetExo = NodeSet(id, convert(Array{Int32}, nsets[name]))

        write_set(exo, nsetExo)
        write_name(exo, nsetExo, name)
    end

    fem_active = !isnothing(topology)
    n_blocks = length(all_block_name_list)

    for (block, block_name) in enumerate(all_block_name_list)
        conn = get_block_nodes(block_Id, block)

        if isempty(conn)
            # An empty block would still consume one of the declared element blocks.
            @warn "Block $block ($block_name) contains no nodes and is skipped. " *
                  "The block count passed to create_result_file must not include it."
            continue
        end

        if fem_active && !isnothing(fem_block) && fem_block[conn[1]]
            # TODO: this writes the complete topology into every FE block. If more than
            # one FE block exists, restrict it to the elements of this block.
            fem_conn = Matrix(topology')
            if dof == 3
                fem_conn[(end - 1):end, :] .= fem_conn[[end; end - 1], :]
                fem_conn[(end - 5):(end - 4), :] .= fem_conn[[end - 4; end - 5], :]
                write_block(exo, block, "HEX8", fem_conn)
            else
                fem_conn[(end - 1):end, :] .= fem_conn[[end; end - 1], :]
                write_block(exo, block, "QUAD4", fem_conn)
            end
            write_name(exo, Block, block, block_name)
        else
            write_block(exo, block, "SPHERE", conn)
            write_name(exo, Block, block, block_name)
        end
    end

    # Bond blocks are written after all material blocks, so their ids never collide
    # with a material block id.
    for (block, bond_block) in bond_blocks
        bond_id = n_blocks + block
        write_block(exo, bond_id, "BAR2", bond_block.conn)
        write_name(exo, Block, bond_id, all_block_name_list[block] * "_bonds")
    end

    # write element id map
    n_bond_elements, _ = bond_export_sizes(bond_blocks)
    if fem_active
        element_ids = Int32.(elem_global_ids)
    else
        element_ids = Int32.(global_ids)
    end
    if n_bond_elements > 0
        # The element map has to cover the bond elements as well, otherwise its length
        # does not match num_elems from the initialisation.
        offset = Int32(maximum(element_ids; init = Int32(0)))
        append!(element_ids, Int32.(offset .+ (1:n_bond_elements)))
    end
    write_id_map(exo, NodeMap, Int32.(global_ids))
    write_id_map(exo, ElementMap, element_ids)

    # output structure var_name -> [fieldname, exodus id, field dof]
    nodal_outputs = Dict(key => value
                         for (key, value) in output["Fields"]
                         if (!value["global_var"] && !get(value, "bond_var", false)))
    global_outputs = Dict(key => value
                          for (key, value) in output["Fields"] if (value["global_var"]))
    nodal_output_names = collect(keys(sort!(OrderedDict(nodal_outputs))))
    global_output_names = collect(keys(sort!(OrderedDict(global_outputs))))

    if length(nodal_output_names) == 0
        @warn "No nodal output variables defined, but exodus file created. Please check your output definition."
    else
        write_names(exo, NodalVariable, nodal_output_names)

        nnodes = length(coords[1, :])
        for varname in nodal_output_names
            write_values(exo, NodalVariable, 1, varname, zeros(Float64, nnodes))
        end
    end

    # Bond variables are element variables on the BAR2 blocks. Exodus.jl has no truth
    # table API, so every element variable has to exist on every element block; the
    # material blocks get zeros.
    if !isempty(bond_output_names)
        if isempty(bond_blocks)
            @warn "Bond output variables are defined, but no bond block was created. " *
                  "Enable \"Bond Export\" for the affected output."
        else
            write_names(exo, ElementVariable, bond_output_names)
            sizes = element_block_sizes(block_Id, n_blocks, bond_blocks)
            for varname in bond_output_names
                for (exo_block_id, n_elements) in sizes
                    write_values(exo, ElementVariable, 1, exo_block_id, varname,
                                 zeros(Float64, n_elements))
                end
            end
        end
    end

    global_used = length(global_output_names) > 0

    if global_used
        write_names(exo, GlobalVariable, global_output_names)
        write_values(exo, GlobalVariable, 1, zeros(Float64, length(global_output_names)))
    end

    return exo
end

"""
    write_step_and_time(exo::ExodusDatabase, step::Int64, time::Float64)

Writes the step and time in the exodus file

# Arguments
- `exo::ExodusDatabase`: The exodus file
- `step::Int64`: The step
- `time::Float64`: The time
# Returns
- `exo::ExodusDatabase`: The exodus file
"""
function write_step_and_time(exo::ExodusDatabase, step::Int64, time::Float64)
    write_time(exo, step, Float64(time))
    return exo
end

"""
    write_nodal_results_in_exodus(exo::ExodusDatabase, step::Int64, output::Dict)

Writes the nodal results in the exodus file

# Arguments
- `exo::ExodusDatabase`: The exodus file
- `step::Int64`: The step
- `output::Dict`: The output
# Returns
- `exo::ExodusDatabase`: The exodus file
"""
function write_nodal_results_in_exodus(exo::ExodusDatabase,
                                       step::Int64,
                                       output::Dict)
    nnodes = Data_Manager.get_nnodes()
    for varname in keys(output)
        field = Data_Manager.get_field(output[varname]["fieldname"],
                                       output[varname]["time"])
        # exo, timestep::Integer, id::Integer, var_index::Integer, vector
        # => https://github.com/cmhamel/Exodus.jl/blob/master/src/Variables.jl
        if haskey(output[varname], "dof")
            var = convert(Array{Float64}, field[1:nnodes, output[varname]["dof"]])
        else
            var = convert(Array{Float64},
                          field[1:nnodes, output[varname]["i_dof"],
                                output[varname]["j_dof"]])
        end
        # interface does not work with Int yet 28//08//2023
        write_values(exo, NodalVariable, step, varname, var)
    end
    return exo
end

"""
    write_bond_results_in_exodus(exo, step, output, bond_blocks, block_Id, n_blocks)

Writes the bond results in the exodus file, as element variables on the BAR2 blocks.

A bond can carry the same variable types as a point (scalar, vector, tensor); the
component selection in `output` therefore follows the same `"dof"` / `"i_dof"` +
`"j_dof"` convention as the nodal output. The values are gathered in the order of the
bond connectivity, so bond `k` of a block always refers to the same element that was
written at initialisation.

Because Exodus.jl offers no truth table, the same variable is also written to the
material blocks — as zeros.

# Arguments
- `exo::ExodusDatabase`: The exodus file
- `step::Int64`: The step
- `output::Dict`: The bond output definitions
- `bond_blocks::AbstractDict{Int64,BondBlock}`: The bond blocks
- `block_Id::AbstractVector{Int64}`: Block id per node
- `n_blocks::Int64`: Number of material blocks
# Returns
- `exo::ExodusDatabase`: The exodus file
"""
function write_bond_results_in_exodus(exo::ExodusDatabase,
                                      step::Int64,
                                      output::Dict,
                                      bond_blocks::AbstractDict{Int64,BondBlock},
                                      block_Id::AbstractVector{Int64},
                                      n_blocks::Int64)
    isempty(output) && return exo
    isempty(bond_blocks) && return exo

    sizes = element_block_sizes(block_Id, n_blocks, bond_blocks)

    for varname in keys(output)
        field = Data_Manager.get_field(output[varname]["fieldname"],
                                       output[varname]["time"])

        if haskey(output[varname], "dof")
            dof = output[varname]["dof"]
        elseif haskey(output[varname], "i_dof")
            dof = (output[varname]["i_dof"], output[varname]["j_dof"])
        else
            dof = nothing
        end

        for (exo_block_id, n_elements) in sizes
            block = exo_block_id - n_blocks
            if haskey(bond_blocks, block)
                values = gather_bond_values(field, bond_blocks[block], dof)
            else
                # material block: the variable is not defined there
                values = zeros(Float64, n_elements)
            end
            write_values(exo, ElementVariable, step, exo_block_id, varname, values)
        end
    end

    return exo
end

"""
    write_global_results_in_exodus(exo::ExodusDatabase, step::Int64, global_values)

Writes the global results in the exodus file

# Arguments
- `exo::ExodusDatabase`: The exodus file
- `step::Int64`: The step
- `global_values`: The global values
# Returns
- `exo::ExodusDatabase`: The exodus file
"""
function write_global_results_in_exodus(exo::ExodusDatabase, step::Int64, global_values)
    write_values(exo, GlobalVariable, step, Vector{Float64}(global_values))
    return exo
end

"""
    merge_exodus_file(file_name::Union{AbstractString,String})

Merges the exodus file

# Arguments
- `file_name::Union{AbstractString,String}`: The name of the file to merge
# Returns
- `exo::ExodusDatabase`: The exodus file
"""
function merge_exodus_file(file_name::Union{AbstractString,String})
    epu(file_name)
end
