# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Material

include("Material_Models/Ordinary/Ordinary.jl")

using TimerOutputs: @timeit
using ....Data_Manager
using ....PeriLabExceptions: @abort
using ....ModuleLoader: find_registered_modules
using .....ParameterSpec: @params, Dependent, register_base!, WithBase, Composite,
                          ParseContext, add_error!, join_path, Table1D, dependent_tables,
                          value, model_parts, model_module,
                          model_function
import .....ParameterSpec: check!

@params "Lowers the yield stress around a point to seed localization (Correspondence Plastic)" struct FlawFunctionParams
    "Switches the flaw on"
    active::Bool = req("Active")
    "Shape of the flaw: a Gaussian dip of the yield stress"
    function_name::String = req("Function"; allowed = ["Pre-defined"])
    "Width of the Gaussian dip around the flaw location"
    flaw_size::Union{Nothing,Float64} = opt("Flaw Size"; default = nothing, min = 0,
                                            quantity = :length)
    "Relative reduction of the yield stress at the flaw location"
    flaw_magnitude::Union{Nothing,Float64} = opt("Flaw Magnitude"; default = nothing,
                                                 min = 0,
                                                 max = 1)
    "x coordinate of the flaw center"
    flaw_location_x::Union{Nothing,Float64} = opt("Flaw Location X"; default = nothing)
    "y coordinate of the flaw center"
    flaw_location_y::Union{Nothing,Float64} = opt("Flaw Location Y"; default = nothing)
    "z coordinate of the flaw center (3D); 0 if not given"
    flaw_location_z::Union{Nothing,Float64} = opt("Flaw Location Z"; default = nothing)
end

# an active flaw function needs its size, magnitude and location (Z is 0 if missing)
function check!(p::FlawFunctionParams, path::String, ctx::ParseContext)
    p.active || return nothing
    required = (p.flaw_size => "Flaw Size", p.flaw_magnitude => "Flaw Magnitude",
                p.flaw_location_x => "Flaw Location X",
                p.flaw_location_y => "Flaw Location Y")
    missing_keys = [alias for (value, alias) in required if value === nothing]
    if !isempty(missing_keys)
        add_error!(ctx, path, "an active Flaw Function needs $(join(missing_keys, ", "))")
        return nothing
    end
    for (value, alias) in required[1:2]
        value > 0 || add_error!(ctx, join_path(path, alias), "must be positive")
    end
    return nothing
end

"""
    MaterialBaseParams

Keys every material model may use (read from the same YAML block as the model's
own keys): symmetry, elastic constants, and options shared by the material
framework. Which elastic constants are needed depends on the model and the
symmetry; they are completed at initialisation.
"""
@params "Material model of the blocks that use this entry: elastic constants, symmetry and options shared by all material models" struct MaterialBaseParams
    "Material symmetry and 2D assumption, e.g. \"isotropic plane strain\", \"isotropic plane stress\", \"orthotropic\", \"transverse isotropic\" or \"anisotropic\"; selects which elastic constants are needed"
    symmetry::Union{Nothing,String} = opt("Symmetry"; default = nothing)
    "Isotropic elastic constant; give two of Young's modulus, Poisson's ratio, bulk modulus and shear modulus, the others are computed"
    youngs_modulus::Union{Nothing,Float64} = opt("Young's Modulus"; default = nothing,
                                                 min = 0,
                                                 quantity = :stress)
    "Isotropic elastic constant; give two of Young's modulus, Poisson's ratio, bulk modulus and shear modulus, the others are computed"
    poissons_ratio::Union{Nothing,Float64} = opt("Poisson's Ratio"; default = nothing,
                                                 min = -1, max = 0.5)
    "Isotropic elastic constant; give two of Young's modulus, Poisson's ratio, bulk modulus and shear modulus, the others are computed"
    bulk_modulus::Union{Nothing,Float64} = opt("Bulk Modulus"; default = nothing, min = 0,
                                               quantity = :stress)
    "Isotropic elastic constant; give two of Young's modulus, Poisson's ratio, bulk modulus and shear modulus, the others are computed"
    shear_modulus::Union{Nothing,Float64} = opt("Shear Modulus"; default = nothing, min = 0,
                                                quantity = :stress)
    "Young's modulus along the material x axis (orthotropic, transverse isotropic); a number or a table over a node field"
    youngs_modulus_x::Union{Nothing,Dependent} = opt("Young's Modulus X"; default = nothing,
                                                     min = 0, quantity = :stress)
    "Young's modulus along the material y axis (orthotropic, transverse isotropic); a number or a table over a node field"
    youngs_modulus_y::Union{Nothing,Dependent} = opt("Young's Modulus Y"; default = nothing,
                                                     min = 0, quantity = :stress)
    "Young's modulus along the material z axis (orthotropic, transverse isotropic); a number or a table over a node field"
    youngs_modulus_z::Union{Nothing,Dependent} = opt("Young's Modulus Z"; default = nothing,
                                                     min = 0, quantity = :stress)
    "Poisson's ratio in the material xy plane (orthotropic, transverse isotropic); a number or a table over a node field"
    poissons_ratio_xy::Union{Nothing,Dependent} = opt("Poisson's Ratio XY";
                                                      default = nothing)
    "Poisson's ratio in the material yz plane (orthotropic, transverse isotropic); a number or a table over a node field"
    poissons_ratio_yz::Union{Nothing,Dependent} = opt("Poisson's Ratio YZ";
                                                      default = nothing)
    "Poisson's ratio in the material xz plane (orthotropic, transverse isotropic); a number or a table over a node field"
    poissons_ratio_xz::Union{Nothing,Dependent} = opt("Poisson's Ratio XZ";
                                                      default = nothing)
    "Shear modulus in the material xy plane (orthotropic, transverse isotropic); a number or a table over a node field"
    shear_modulus_xy::Union{Nothing,Dependent} = opt("Shear Modulus XY"; default = nothing,
                                                     min = 0, quantity = :stress)
    "Shear modulus in the material yz plane (orthotropic, transverse isotropic); a number or a table over a node field"
    shear_modulus_yz::Union{Nothing,Dependent} = opt("Shear Modulus YZ"; default = nothing,
                                                     min = 0, quantity = :stress)
    "Shear modulus in the material xz plane (orthotropic, transverse isotropic); a number or a table over a node field"
    shear_modulus_xz::Union{Nothing,Dependent} = opt("Shear Modulus XZ"; default = nothing,
                                                     min = 0, quantity = :stress)
    "Entry 11 of the stiffness matrix in Voigt notation (anisotropic)"
    c11::Union{Nothing,Float64} = opt("C11"; default = nothing, quantity = :stress)
    "Entry 12 of the stiffness matrix in Voigt notation (anisotropic)"
    c12::Union{Nothing,Float64} = opt("C12"; default = nothing, quantity = :stress)
    "Entry 13 of the stiffness matrix in Voigt notation (anisotropic)"
    c13::Union{Nothing,Float64} = opt("C13"; default = nothing, quantity = :stress)
    "Entry 14 of the stiffness matrix in Voigt notation (anisotropic)"
    c14::Union{Nothing,Float64} = opt("C14"; default = nothing, quantity = :stress)
    "Entry 15 of the stiffness matrix in Voigt notation (anisotropic)"
    c15::Union{Nothing,Float64} = opt("C15"; default = nothing, quantity = :stress)
    "Entry 16 of the stiffness matrix in Voigt notation (anisotropic)"
    c16::Union{Nothing,Float64} = opt("C16"; default = nothing, quantity = :stress)
    "Entry 22 of the stiffness matrix in Voigt notation (anisotropic)"
    c22::Union{Nothing,Float64} = opt("C22"; default = nothing, quantity = :stress)
    "Entry 23 of the stiffness matrix in Voigt notation (anisotropic)"
    c23::Union{Nothing,Float64} = opt("C23"; default = nothing, quantity = :stress)
    "Entry 24 of the stiffness matrix in Voigt notation (anisotropic)"
    c24::Union{Nothing,Float64} = opt("C24"; default = nothing, quantity = :stress)
    "Entry 25 of the stiffness matrix in Voigt notation (anisotropic)"
    c25::Union{Nothing,Float64} = opt("C25"; default = nothing, quantity = :stress)
    "Entry 26 of the stiffness matrix in Voigt notation (anisotropic)"
    c26::Union{Nothing,Float64} = opt("C26"; default = nothing, quantity = :stress)
    "Entry 33 of the stiffness matrix in Voigt notation (anisotropic)"
    c33::Union{Nothing,Float64} = opt("C33"; default = nothing, quantity = :stress)
    "Entry 34 of the stiffness matrix in Voigt notation (anisotropic)"
    c34::Union{Nothing,Float64} = opt("C34"; default = nothing, quantity = :stress)
    "Entry 35 of the stiffness matrix in Voigt notation (anisotropic)"
    c35::Union{Nothing,Float64} = opt("C35"; default = nothing, quantity = :stress)
    "Entry 36 of the stiffness matrix in Voigt notation (anisotropic)"
    c36::Union{Nothing,Float64} = opt("C36"; default = nothing, quantity = :stress)
    "Entry 44 of the stiffness matrix in Voigt notation (anisotropic)"
    c44::Union{Nothing,Float64} = opt("C44"; default = nothing, quantity = :stress)
    "Entry 45 of the stiffness matrix in Voigt notation (anisotropic)"
    c45::Union{Nothing,Float64} = opt("C45"; default = nothing, quantity = :stress)
    "Entry 46 of the stiffness matrix in Voigt notation (anisotropic)"
    c46::Union{Nothing,Float64} = opt("C46"; default = nothing, quantity = :stress)
    "Entry 55 of the stiffness matrix in Voigt notation (anisotropic)"
    c55::Union{Nothing,Float64} = opt("C55"; default = nothing, quantity = :stress)
    "Entry 56 of the stiffness matrix in Voigt notation (anisotropic)"
    c56::Union{Nothing,Float64} = opt("C56"; default = nothing, quantity = :stress)
    "Entry 66 of the stiffness matrix in Voigt notation (anisotropic)"
    c66::Union{Nothing,Float64} = opt("C66"; default = nothing, quantity = :stress)
    "Column of the state variables that scales the isotropic elastic constants node by node, e.g. for a degradation computed in a UMAT"
    state_factor_id::Union{Nothing,Int64} = opt("State Factor ID"; default = nothing,
                                                min = 1)
    "Polynomial order of the bond-associated correspondence formulation"
    accuracy_order::Union{Nothing,Int64} = opt("Accuracy Order"; default = nothing, min = 1)
    "Stabilization against zero-energy modes of correspondence materials, e.g. \"Global\""
    zero_energy_control::Union{Nothing,String} = opt("Zero Energy Control";
                                                     default = nothing)
    "Uses the bond-associated correspondence formulation, which evaluates the deformation gradient per bond instead of per node"
    bond_associated::Bool = opt("Bond Associated"; default = false)
    "Uses the small-strain tensor instead of the Green-Lagrange strain in correspondence materials"
    linear_strain::Bool = opt("Linear Strain"; default = false)
    flaw_function::Union{Nothing,FlawFunctionParams} = opt("Flaw Function";
                                                           default = nothing)
end

const _ORTHOTROPIC_KEYS = ((:youngs_modulus_x, "Young's Modulus X"),
                           (:youngs_modulus_y, "Young's Modulus Y"),
                           (:youngs_modulus_z, "Young's Modulus Z"),
                           (:poissons_ratio_xy, "Poisson's Ratio XY"),
                           (:poissons_ratio_yz, "Poisson's Ratio YZ"),
                           (:poissons_ratio_xz, "Poisson's Ratio XZ"),
                           (:shear_modulus_xy, "Shear Modulus XY"),
                           (:shear_modulus_yz, "Shear Modulus YZ"),
                           (:shear_modulus_xz, "Shear Modulus XZ"))

function _transverse_keys(symmetry::String)
    keys = [(:youngs_modulus_x, "Young's Modulus X"),
        (:youngs_modulus_y, "Young's Modulus Y"),
        (:poissons_ratio_xy, "Poisson's Ratio XY")]
    occursin("plane stress", symmetry) ||
        push!(keys, (:poissons_ratio_yz, "Poisson's Ratio YZ"))
    push!(keys, (:shear_modulus_xy, "Shear Modulus XY"))
    if !occursin("plane strain", symmetry) && !occursin("plane stress", symmetry)
        push!(keys, (:shear_modulus_yz, "Shear Modulus YZ"))
    end
    return keys
end

# stiffness-matrix materials must define all their constants
function check!(p::MaterialBaseParams, path::String, ctx::ParseContext)
    p.symmetry === nothing && return nothing
    symmetry = lowercase(p.symmetry)
    required = if occursin("anisotropic", symmetry)
        [(Symbol("c$i$j"), "C$i$j") for i in 1:6 for j in i:6]
    elseif occursin("transverse isotropic", symmetry)
        _transverse_keys(symmetry)
    elseif occursin("orthotropic", symmetry)
        collect(_ORTHOTROPIC_KEYS)
    else
        Tuple{Symbol,String}[]
    end
    missing_keys = [alias for (field, alias) in required if getfield(p, field) === nothing]
    isempty(missing_keys) ||
        add_error!(ctx, join_path(path, "Symmetry"),
                   "\"$(p.symmetry)\" requires $(join(missing_keys, ", "))")
    return nothing
end

# registration runs at load time, never during precompilation
__init__() = register_base!(:material, MaterialBaseParams)

# the correspondence models register as materials too; they live in Correspondence
const CORRESPONDENCE_DIR = joinpath(@__DIR__, "Material_Models", "Correspondence")
const CORRESPONDENCE_TEMPLATE = joinpath(@__DIR__, "Material_Models", "Material_template",
                                         "correspondence_template.jl")
for file in find_registered_modules(@__DIR__, "register_material";
                                    exclude = [CORRESPONDENCE_DIR, CORRESPONDENCE_TEMPLATE])
    include(file)
end
include("Material_Models/Correspondence/Correspondence.jl")

using ...Material_Basis:
                         distribute_forces!,
                         init_local_damping_due_to_damage,
                         local_damping_due_to_damage
using LinearAlgebra: dot
using StaticArrays
"""
    ElasticModuli

Completed isotropic elastic constants of a block. Each entry is a `Float64`, or
a node field (`Vector{Float64}`) when moduli are given per node (mesh columns
`Bulk_Modulus`, …) or scaled by a state variable. Read one with `modulus(x, iID)`.
"""
struct ElasticModuli{K,E,G,N}
    bulk_modulus::K
    youngs_modulus::E
    shear_modulus::G
    poissons_ratio::N
end

@inline modulus(x::Real, ::Int64) = x
@inline modulus(x::AbstractVector, iID::Int64) = x[iID]

"""
    BlockMaterial

Typed material of one block: the shared base parameters, the model struct (or
`Composite`), the symmetry used by the force models (`"plane strain"`,
`"plane stress"` or `"3D"`), and the completed elastic moduli (`nothing` for
materials defined by a stiffness matrix).
"""
# mutable only for heap identity (never mutated): passing it through dynamic calls
# then shares the object instead of copying the struct
mutable struct BlockMaterial{B,M,E}
    base::B
    model::M
    symmetry::String          # used by the force models: "plane strain", "plane stress" or "3D"
    hooke_symmetry::String    # used by the Hooke matrix (see hooke_symmetry)
    moduli::E
    tables::Vector{Table1D}   # dependent tables of base and model, re-bound every step
    extras::Dict{String,Any}  # values of indexed keys, e.g. Property_1
    correspondence::Bool      # the model parts are correspondence models
end

"""
    material_symmetry(symmetry, dof)

The symmetry the force models use. Plane strain / plane stress are ignored in 3D;
every other symmetry gives `"3D"`.
"""
function material_symmetry(symmetry::Union{Nothing,String}, dof::Int64)
    symmetry === nothing && return "3D"
    s = symmetry
    if dof == 3
        s = replace(replace(s, r"plane strain$" => ""), r"plane stress$" => "")
    end
    s = lowercase(s)
    occursin("plane strain", s) && return "plane strain"
    occursin("plane stress", s) && return "plane stress"
    return "3D"
end

"""
    hooke_symmetry(symmetry, dof)

The symmetry string the Hooke matrix uses: as given, with a trailing plane strain /
plane stress removed in 3D, `"isotropic"` if missing.
"""
function hooke_symmetry(symmetry::Union{Nothing,String}, dof::Int64)
    symmetry === nothing && return "isotropic"
    dof == 3 || return symmetry
    return replace(replace(symmetry, r"plane strain$" => ""), r"plane stress$" => "")
end

const _MODULI = (("Bulk Modulus", :bulk_modulus), ("Young's Modulus", :youngs_modulus),
                 ("Shear Modulus", :shear_modulus), ("Poisson's Ratio", :poissons_ratio))

_modulus_field(key::String) = replace(key, " " => "_")

# value of one modulus: the node field if it exists, a new node field if any
# modulus is per node, otherwise the given constant (0.0 if not given)
function _modulus_value(given, field_allocated::Bool, any_field_allocated::Bool,
                        key::String)
    field_allocated && return Data_Manager.get_field(_modulus_field(key))
    if any_field_allocated
        return given === nothing ?
               Data_Manager.create_constant_node_scalar_field(_modulus_field(key),
                                                              Float64) :
               Data_Manager.create_constant_node_scalar_field(_modulus_field(key), Float64;
                                                              default_value = given)
    end
    return given === nothing ? 0.0 : given
end

"""
    elastic_moduli(base, bond_based, dof)

Completes the isotropic elastic constants from any two of bulk modulus, Young's
modulus, shear modulus and Poisson's ratio (bond-based models: Poisson's ratio
is fixed). Moduli given per node (fields `Bulk_Modulus`, …) or a `State Factor
ID` make the result per node, and the node fields are updated. Returns `nothing`
for anisotropic, orthotropic and transverse isotropic materials (stiffness
matrix; completeness is checked when the input is read).
"""
function elastic_moduli(base::MaterialBaseParams, bond_based::Bool, dof::Int64)
    state_factor_defined = base.state_factor_id !== nothing
    allocated = Dict(key => Data_Manager.has_key(_modulus_field(key))
                     for (key, _) in _MODULI)
    any_field_allocated = any(values(allocated)) || state_factor_defined
    given = Dict(key => getfield(base, name) for (key, name) in _MODULI)
    has = Dict(key => given[key] !== nothing || allocated[key] for (key, _) in _MODULI)

    K = _modulus_value(given["Bulk Modulus"], allocated["Bulk Modulus"],
                       any_field_allocated, "Bulk Modulus")
    E = _modulus_value(given["Young's Modulus"], allocated["Young's Modulus"],
                       any_field_allocated, "Young's Modulus")
    G = _modulus_value(given["Shear Modulus"], allocated["Shear Modulus"],
                       any_field_allocated, "Shear Modulus")
    nu = _modulus_value(given["Poisson's Ratio"], allocated["Poisson's Ratio"],
                        any_field_allocated, "Poisson's Ratio")
    bulk, youngs, shear, poissons = has["Bulk Modulus"], has["Young's Modulus"],
                                    has["Shear Modulus"], has["Poisson's Ratio"]

    if bond_based
        nu_fixed = dof == 2 ? 1 / 3 : 1 / 4
        if nu != 0.0 && nu != nu_fixed
            @warn "Chosen Bond-based model only supports a fixed Poisson's ratio of " *
                  string(nu_fixed)
        end
        nu = nu_fixed
        poissons = true
    end
    if base.symmetry !== nothing
        symmetry = lowercase(base.symmetry)
        if occursin("anisotropic", symmetry) ||
           occursin("transverse isotropic", symmetry) ||
           occursin("orthotropic", symmetry)
            return nothing
        end
    else
        @warn "Material symmetry is not defined, assuming isotropic material"
    end

    if bulk + youngs + shear + poissons < 2
        @abort "Minimum of two parameters are needed for isotropic material"
    elseif bulk + youngs + shear + poissons > 2
        @warn "Only two parameters are needed for isotropic material, ignoring additional parameters"
    end

    if bulk && poissons
        E = 3 .* K .* (1 .- 2 .* nu)
        G = 3 .* K .* (1 .- 2 .* nu) ./ (2 .+ 2 .* nu)
    end
    if shear && poissons
        E = 2 .* G .* (1 .+ nu)
        K = 2 .* G .* (1 .+ nu) ./ (3 .- 6 .* nu)
    end
    if bulk && shear
        E = 9 .* K .* G ./ (3 .* K .+ G)
        nu = (3 .* K .- 2 .* G) ./ (6 .* K .+ 2 .* G)
    end
    if youngs && shear
        K = E .* G ./ (9 .* G .- 3 .* E)
        nu = E ./ (2 .* G) .- 1
    end
    if youngs && bulk
        G = 3 .* K .* E ./ (9 .* K .- E)
        nu = (3 .* K .- E) ./ (6 .* K)
    end
    if youngs && poissons
        K = E ./ (3 .- 6 .* nu)
        G = E ./ (2 .+ 2 .* nu)
    end

    if state_factor_defined && Data_Manager.has_key("State Variables")
        state_factor = Data_Manager.get_field("State Variables")[:, base.state_factor_id]
        K = K .* state_factor
        E = E .* state_factor
        G = G .* state_factor
    end
    if any_field_allocated
        Data_Manager.get_field("Bulk_Modulus") .= K
        Data_Manager.get_field("Young's_Modulus") .= E
        Data_Manager.get_field("Shear_Modulus") .= G
        Data_Manager.get_field("Poisson's_Ratio") .= nu
    end
    return ElasticModuli(K, E, G, nu)
end

"""
    block_material(wb, dof)

The `BlockMaterial` of a parsed material block (bond-based models fix Poisson's
ratio).
"""
function block_material(wb::WithBase, dof::Int64)
    model_name = wb.name
    bond_based = occursin("Bond-based", model_name) &&
                 !occursin("Unified Bond-based", model_name)
    correspondence = correspondence_model(wb.model, model_name)
    moduli = elastic_moduli(wb.base, bond_based, dof)
    if moduli === nothing && !correspondence
        @abort "Material Model \"$model_name\": Symmetry \"$(wb.base.symmetry)\" defines a stiffness matrix, which only correspondence models support."
    end
    return BlockMaterial(wb.base, wb.model, material_symmetry(wb.base.symmetry, dof),
                         hooke_symmetry(wb.base.symmetry, dof), moduli,
                         vcat(dependent_tables(wb.base), dependent_tables(wb.model)),
                         wb.extras,
                         correspondence)
end

"A model struct defined in the Correspondence module (a correspondence formulation)."
is_correspondence(part) = parentmodule(model_module(part)) === Correspondence

"""
    correspondence_model(model, name)

`true` if the parts of `model` are correspondence models. Correspondence and other
material models cannot be combined with `+`.
"""
function correspondence_model(model, name::AbstractString)
    parts = model_parts(model)
    correspondence = any(is_correspondence, parts)
    if correspondence && !all(is_correspondence, parts)
        @abort "Material Model \"$name\": correspondence models cannot be combined with other material models."
    end
    return correspondence
end

_constant(x) = x === nothing ? nothing : value(x, 1)

"""
    critical_bulk_modulus(material)

Bulk modulus used for the critical time step (legacy rules): the isotropic bulk
modulus; for orthotropic constants the compliance-based estimate; else C44/C55/C66
or Shear Modulus XY as estimates; `nothing` if none is defined.
"""
function critical_bulk_modulus(material::BlockMaterial)
    material.moduli === nothing || return material.moduli.bulk_modulus
    base = material.base
    nu_xy, nu_yz, nu_xz = _constant(base.poissons_ratio_xy),
                          _constant(base.poissons_ratio_yz),
                          _constant(base.poissons_ratio_xz)
    if nu_xy !== nothing && nu_yz !== nothing && nu_xz !== nothing
        E_x, E_y, E_z = _constant(base.youngs_modulus_x), _constant(base.youngs_modulus_y),
                        _constant(base.youngs_modulus_z)
        s11 = 1 / E_x
        s22 = 1 / E_y
        s33 = 1 / E_z
        s12 = -nu_xy / E_x
        s23 = -nu_yz / E_z
        s13 = -nu_xz / E_z
        return 1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s13))
    elseif base.c44 !== nothing && base.c55 !== nothing && base.c66 !== nothing
        return maximum([base.c44 / 2, base.c55 / 2, base.c66 / 2])
    elseif base.shear_modulus_xy !== nothing
        return _constant(base.shear_modulus_xy) / 2
    end
    return nothing
end

export init_model
export compute_model
export distribute_force_densities
export init_fields
export fields_for_local_synchronization
export compute_local_damping
export init_local_damping

function compute_local_damping(nodes, params, dt)
    return local_damping_due_to_damage(nodes, params, dt)
end
function init_local_damping(nodes, symmetry::String, local_damping)
    return init_local_damping_due_to_damage(nodes, symmetry, local_damping)
end

"""
    init_fields()

Initialize material model fields
"""
function init_fields()
    dof = Data_Manager.get_dof()
    Data_Manager.create_node_vector_field("Forces", Float64, dof) #-> only if it is an output
    # tbd later in the compute class
    Data_Manager.create_constant_node_vector_field("External Forces", Float64, dof)
    Data_Manager.create_node_vector_field("Force Densities", Float64, dof)
    Data_Manager.create_constant_node_vector_field("External Force Densities", Float64, dof)
    Data_Manager.create_node_vector_field("Acceleration", Float64, dof)
    Data_Manager.create_node_vector_field("Velocity", Float64, dof)
    Data_Manager.create_constant_bond_vector_state("Bond Forces", Float64, dof)
    Data_Manager.create_constant_bond_scalar_state("Temporary Bond Field", Float64)
    Data_Manager.create_node_vector_field("Displacements", Float64, dof)
    Data_Manager.create_bond_vector_state("Deformed Bond Geometry", Float64, dof)
    Data_Manager.create_bond_scalar_state("Deformed Bond Length", Float64)
    # Data_Manager.set_synch("Bond Forces", false, true)
    Data_Manager.set_synch("Force Densities", true, false)
    Data_Manager.set_synch("Velocity", false, true)
    Data_Manager.set_synch("Displacements", false, true)
    Data_Manager.set_synch("Acceleration", false, true)
    Data_Manager.set_synch("Deformed Coordinates", false, true)
end

"""
    init_model(nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the material model.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    material = Data_Manager.get_block_models(block).material
    if material === nothing
        @abort "Block " * string(block) * " has no material model defined."
        return
    end

    bind_material!(material)
    if material.correspondence
        return Correspondence.init_model(nodes, material, block)
    end

    for part in model_parts(material.model)
        model_module(part).init_model(nodes, part, material, block)
    end
    #TODO in extra function
    # nlist = Data_Manager.get_nlist()
    nlist_filtered_ids = Data_Manager.get_filtered_nlist()
    if !isnothing(nlist_filtered_ids)
        bond_norm = Data_Manager.get_field("Bond Norm")
        bond_geometry = Data_Manager.get_field("Bond Geometry")
        for iID in nodes
            if length(nlist_filtered_ids[iID]) != 0
                for neighborID in nlist_filtered_ids[iID]
                    bond_norm[iID][neighborID] = bond_norm[iID][neighborID] *
                                                 sign(dot((bond_geometry[iID][neighborID]),
                                                          bond_norm[iID][neighborID]))
                end
            end
        end
    end
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block id
"""
function fields_for_local_synchronization(model, block)
    material = Data_Manager.get_block_models(block).material
    if material.correspondence
        return Correspondence.fields_for_local_synchronization(model, material)
    end

    for part in model_parts(material.model)
        model_module(part).fields_for_local_synchronization(model)
    end
end

"""
    compute_model(nodes::AbstractVector{Int64}, material, block::Int64, time::Float64, dt::Float64)

Computes the material models

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `material::BlockMaterial`: The typed block material
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64},
                       material,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    @timeit "all" begin
        if material.correspondence
            @timeit "corresponcence" begin
                Correspondence.compute_model(nodes, bind_material!(material), block, time,
                                             dt)
                return
            end
        end
        @timeit "material" compute_block_material(nodes, material, block, time, dt)
    end
end

"""
    bind_material!(material)

Binds the dependent tables of a block material to the current node fields. Call
it before every evaluation, because the N/NP1 field arrays are swapped every step.
"""
function bind_material!(material::BlockMaterial)
    Data_Manager.bind_dependent_tables!(material.tables)
    return material
end

# function barrier: `material` has a concrete type here
function compute_block_material(nodes::AbstractVector{Int64}, material::BlockMaterial,
                                block::Int64, time::Float64, dt::Float64)
    bind_material!(material)
    for part in model_parts(material.model)
        model_function(part, Val(:compute_model))(nodes, part, material, block, time, dt)
    end
    return nothing
end

"""
    check_material_symmetry(material, dof)

2D needs plane strain or plane stress in `Symmetry`; in 3D these are ignored
(with a warning). A missing Symmetry is not checked here (isotropic).
"""
function check_material_symmetry(material::BlockMaterial, dof::Int64)
    symmetry = material.base.symmetry
    symmetry === nothing && return nothing
    if dof == 2 && !occursin("plane strain", symmetry) &&
       !occursin("plane stress", symmetry)
        @abort "Model definition is missing; plane stress or plane strain has to be defined for 2D"
        return
    end
    if dof == 3 && occursin("plane strain", symmetry)
        @warn "Plane strain symmetry is not supported for 3D, going to ignore it"
    end
    if dof == 3 && occursin("plane stress", symmetry)
        @warn "Plane stress symmetry is not supported for 3D, going to ignore it"
    end
    return nothing
end

"""
    distribute_force_densities(nodes::AbstractVector{Int64})

Distribute the force densities.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
"""
function distribute_force_densities(nodes::AbstractVector{Int64})
    @timeit "load data" begin
        nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
        nlist_filtered_ids = Data_Manager.get_filtered_nlist()
        bond_force::BondVectorField{Float64} = Data_Manager.get_field("Bond Forces")
        force_densities::NodeVectorField{Float64} = Data_Manager.get_field("Force Densities",
                                                                           "NP1")
        volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
        bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    end
    if !isnothing(nlist_filtered_ids)
        bond_norm = Data_Manager.get_field("Bond Norm")
        displacements = Data_Manager.get_field("Displacements", "NP1")
        @timeit "local dist" distribute_forces!(force_densities,
                                                nodes,
                                                nlist,
                                                nlist_filtered_ids,
                                                bond_force,
                                                volume,
                                                bond_damage,
                                                displacements,
                                                bond_norm)
    else
        @timeit "local dist" distribute_forces!(force_densities, nodes, nlist,
                                                bond_force, volume, bond_damage)
    end
end

function compute_correspondence_bond_forces(nodes::AbstractVector{Int64}, material,
                                            block::Int64, time::Float64, dt::Float64)
    Correspondence.compute_correspondence_model(nodes, bind_material!(material), block,
                                                time, dt)
end

end
