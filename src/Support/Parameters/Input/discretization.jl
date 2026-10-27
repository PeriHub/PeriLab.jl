# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params "Finite element topology read from a separate file, for FEM or coupled FEM/PD meshes" struct ExternalTopologyParams
    "File with one element per line, listing the node numbers of the element"
    file::String = req("File")
    "Keeps the peridynamic neighbors of nodes in finite elements in addition to the element neighbors"
    add_neighbor_search::Union{Nothing,Bool} = opt("Add Neighbor Search"; default = nothing)
end

@params "Adds two blocks of regularly spaced nodes on both sides of the mesh, e.g. as load introduction" struct SurfaceExtrusionParams
    "Axis along which the mesh is extended on both sides"
    direction::String = req("Direction"; allowed = ["X", "Y", "Z"])
    "Node spacing along the extrusion direction"
    step_x::Float64 = req("Step_X"; quantity = :length)
    "Node spacing across the extrusion direction in the x-y plane"
    step_y::Float64 = req("Step_Y"; quantity = :length)
    "Node spacing of the added nodes in z (3D only)"
    step_z::Float64 = req("Step_Z"; quantity = :length)
    "Number of node layers added on each side"
    number::Int64 = req("Number"; min = 0)
end

@params "Generates the mesh from a G-code file by sampling the printed paths" struct GcodeParams
    "Regenerate the mesh from the G-code; if false and the .txt mesh written by a previous run exists next to the G-code file, that mesh is read instead"
    overwrite_mesh::Bool = opt("Overwrite Mesh"; default = true)
    "Distance between two nodes along a printed path"
    sampling::Float64 = req("Sampling"; min = 0, quantity = :length)
    "Width of the printed strand; with sampling and height it sets the node volume"
    width::Float64 = req("Width"; min = 0, quantity = :length)
    "Layer height of the printed strand; with sampling and width it sets the node volume"
    height::Float64 = req("Height"; min = 0, quantity = :length)
    "Comma-separated G-code commands that switch node generation on; without it every extruding move creates nodes"
    start_command::Union{Nothing,String} = opt("Start Command"; default = nothing)
    "Comma-separated G-code commands that switch node generation off; needs a start command"
    stop_command::Union{Nothing,String} = opt("Stop Command"; default = nothing)
    "G-code command that ends reading the file"
    end_command::Union{Nothing,String} = opt("End Command"; default = nothing)
    "Block id => condition on the node position x, y, z (e.g. \"z > 2.0\"); a node gets the last block whose condition holds, block 1 otherwise"
    blocks::Union{Nothing,Dict{String,String}} = opt("Blocks"; default = nothing,
                                                     evaluated = true)
end

const _DISK = Dict("Type" => ["Disk"])
const _RECTANGLE = Dict("Type" => ["Rectangular_Plane"])

@params "Removes the bonds crossing a disk or rectangle, e.g. to model a pre-crack" struct BondFilterParams
    "Shape of the filter"
    type::String = req("Type")
    "x component of the normal of the filter plane"
    normal_x::Float64 = req("Normal X")
    "y component of the normal of the filter plane"
    normal_y::Float64 = req("Normal Y")
    "z component of the normal of the filter plane (3D only)"
    normal_z::Union{Nothing,Float64} = opt("Normal Z"; default = nothing)
    "x coordinate of the disk center"
    center_x::Union{Nothing,Float64} = opt("Center X"; default = nothing, when = _DISK)
    "y coordinate of the disk center"
    center_y::Union{Nothing,Float64} = opt("Center Y"; default = nothing, when = _DISK)
    "z coordinate of the disk center (3D only)"
    center_z::Union{Nothing,Float64} = opt("Center Z"; default = nothing, when = _DISK)
    "Radius of the disk"
    radius::Union{Nothing,Float64} = opt("Radius"; default = nothing, min = 0, when = _DISK)
    "x coordinate of the lower left corner of the rectangle"
    lower_left_corner_x::Union{Nothing,Float64} = opt("Lower Left Corner X";
                                                      default = nothing,
                                                      when = _RECTANGLE)
    "y coordinate of the lower left corner of the rectangle"
    lower_left_corner_y::Union{Nothing,Float64} = opt("Lower Left Corner Y";
                                                      default = nothing,
                                                      when = _RECTANGLE)
    "z coordinate of the lower left corner of the rectangle (3D only)"
    lower_left_corner_z::Union{Nothing,Float64} = opt("Lower Left Corner Z";
                                                      default = nothing,
                                                      when = _RECTANGLE)
    "x component of the direction of the rectangle's bottom edge"
    bottom_unit_vector_x::Union{Nothing,Float64} = opt("Bottom Unit Vector X";
                                                       default = nothing, when = _RECTANGLE)
    "y component of the direction of the rectangle's bottom edge"
    bottom_unit_vector_y::Union{Nothing,Float64} = opt("Bottom Unit Vector Y";
                                                       default = nothing, when = _RECTANGLE)
    "z component of the direction of the rectangle's bottom edge (3D only)"
    bottom_unit_vector_z::Union{Nothing,Float64} = opt("Bottom Unit Vector Z";
                                                       default = nothing, when = _RECTANGLE)
    "Length of the rectangle along the bottom unit vector"
    bottom_length::Union{Nothing,Float64} = opt("Bottom Length"; default = nothing, min = 0,
                                                when = _RECTANGLE)
    "Length of the rectangle perpendicular to the bottom edge"
    side_length::Union{Nothing,Float64} = opt("Side Length"; default = nothing, min = 0,
                                              when = _RECTANGLE)
    "Keeps the filtered bonds as contact bonds so the two sides cannot penetrate each other"
    allow_contact::Bool = opt("Allow Contact"; default = false)
end

# the mesh importers (`mesh_import_name()` of the modules in IO/Mesh_Import); an unknown
# type is reported when the mesh is read
const _MESH_TYPES = ["Text File", "Exodus", "Abaqus", "Gmsh", "Gcode"]

@params "Mesh input: the point cloud or mesh file, node sets, bond filters and the neighbor search" struct DiscretizationParams
    "Mesh format"
    type::String = req("Type")
    "Mesh file, relative to the input deck"
    input_mesh_file::String = req("Input Mesh File")
    "Named node sets: a node id, ids separated by spaces, a range \"1:10\", a .txt file with node ids, a condition on x, y, z (e.g. \"x < 0.1\") or \"All\"; Exodus meshes bring their own node sets"
    node_sets::Dict{String,
                    Union{Int64,String}} = opt("Node Sets";
                                               default = Dict{String,
                                                              Union{Int64,String}}(),
                                               evaluated = true,
                                               when = Dict("Type" => ["Text File", "Gmsh",
                                                               "Gcode"]))
    "Named bond filters"
    bond_filters::Dict{String,BondFilterParams} = opt("Bond Filters";
                                                      default = Dict{String,
                                                                     BondFilterParams}())
    gcode::Union{Nothing,GcodeParams} = opt("Gcode"; default = nothing,
                                            when = Dict("Type" => ["Gcode"]))
    "Weight of each bond as a function of its length xi (or xiX, xiY, xiZ), e.g. \"1/xi^2\"; all bonds weigh 1 if not given"
    influence_function::Union{Nothing,String} = opt("Influence Function"; default = nothing,
                                                    evaluated = true)
    surface_extrusion::Union{Nothing,SurfaceExtrusionParams} = opt("Surface Extrusion";
                                                                   default = nothing)
    input_external_topology::Union{Nothing,
                                   ExternalTopologyParams} = opt("Input External Topology";
                                                                 default = nothing)
    "Factor on the horizon in x for the neighbor search, for elliptic neighborhoods"
    horizon_mesh_scaling_x::Union{Nothing,Float64} = opt("Horizon Mesh Scaling X";
                                                         default = nothing)
    "Factor on the horizon in y for the neighbor search, for elliptic neighborhoods"
    horizon_mesh_scaling_y::Union{Nothing,Float64} = opt("Horizon Mesh Scaling Y";
                                                         default = nothing)
    "Factor on the horizon in z for the neighbor search, for elliptic neighborhoods"
    horizon_mesh_scaling_z::Union{Nothing,Float64} = opt("Horizon Mesh Scaling Z";
                                                         default = nothing)
    "How nodes are split across MPI ranks: \"Neighbor based\", \"Node based\", or any other value for an even split by node number"
    distribution_type::Union{Nothing,String} = opt("Distribution Type"; default = nothing)
end

"Horizon scaling per direction, 1.0 where not given."
function mesh_scaling(d::DiscretizationParams)
    return [something(d.horizon_mesh_scaling_x, 1.0),
        something(d.horizon_mesh_scaling_y, 1.0),
        something(d.horizon_mesh_scaling_z, 1.0)]
end

function check!(g::GcodeParams, path::String, ctx::ParseContext)
    g.blocks === nothing && return nothing
    for key in keys(g.blocks)
        tryparse(Int64, key) === nothing &&
            add_error!(ctx, join_path(join_path(path, "Blocks"), key),
                       "expected a block id (integer) as key")
    end
    return nothing
end

"Gcode block assignment: block id => condition, or `nothing`."
function gcode_block_ids(g::GcodeParams)
    g.blocks === nothing && return nothing
    return Dict{Int64,String}(parse(Int64, key) => condition
                              for (key, condition) in g.blocks)
end

"A registered bond filter: the YAML keys its Type needs and the function that applies it."
struct BondFilter
    required::Vector{String}
    run::Function
end

# filter type => BondFilter; only mutated at runtime (module `__init__` functions)
const BOND_FILTERS = Dict{String,BondFilter}()

"""
    register_bond_filter(name, run; required = String[])

Makes bond filter `name` available as `Type` of a `Bond Filters` entry. `run`
is called as `run(nnodes, data, filter::BondFilterParams, nlist, dof)` and
returns the bond flags and the filter normal. `required` lists the
`BondFilterParams` keys the filter needs (Z components only in 3D are checked
by the filter itself). Call it from the filter module's `__init__()`.
"""
function register_bond_filter(name::AbstractString, run::Function;
                              required = String[])
    keys = Set(spec.alias for spec in ParameterSpec.parameter_spec(BondFilterParams))
    for key in required
        key in keys ||
            throw(ParameterSpec.ParamsDefinitionError("bond filter \"$name\": \"$key\" is not a Bond Filters key"))
    end
    BOND_FILTERS[String(name)] = BondFilter(collect(String, required), run)
    return nothing
end

ParameterSpec.schema_choices(::Type{DiscretizationParams}) = Dict("Type" => _MESH_TYPES)
function ParameterSpec.schema_choices(::Type{BondFilterParams})
    Dict("Type" => sort!(collect(keys(BOND_FILTERS))))
end

function ParameterSpec.conditional_required(::Type{BondFilterParams})
    return [Dict("Type" => Any[name]) => BOND_FILTERS[name].required
            for name in sort!(collect(keys(BOND_FILTERS)))
            if !isempty(BOND_FILTERS[name].required)]
end

"The registered bond filter `name`, or `nothing`."
bond_filter(name::AbstractString) = get(BOND_FILTERS, name, nothing)

function check!(f::BondFilterParams, path::String, ctx::ParseContext)
    filter = bond_filter(f.type)
    if filter === nothing
        available = sort!(collect(keys(BOND_FILTERS)))
        suggestion = ParameterSpec.suggest(f.type, available)
        add_error!(ctx, join_path(path, "Type"),
                   suggestion === nothing ?
                   "bond filter \"$(f.type)\" not found; available: $(join(available, ", "))" :
                   "bond filter \"$(f.type)\" not found — did you mean \"$suggestion\"?")
        return nothing
    end
    field = Dict(spec.alias => spec.name
                 for spec in ParameterSpec.parameter_spec(BondFilterParams))
    missing = [key for key in filter.required if getfield(f, field[key]) === nothing]
    isempty(missing) ||
        add_error!(ctx, path, "\"$(f.type)\" bond filter requires: $(join(missing, ", "))")
    return nothing
end
