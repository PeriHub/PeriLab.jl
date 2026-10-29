# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Test-only copies of the legacy (material dict) functions removed when the typed input replaced the material dicts.
# The typed implementation is compared against them, so its numbers stay pinned.
module LegacyMaterialOracle
using PeriLab.Data_Manager
using PeriLab.PeriLabExceptions: @abort
using PeriLab.ParameterSpec: evaluate
using PeriLab.Solver_Manager.Material_Basis: _hooke_matrix

# verbatim copies of the legacy dependent-value helpers (deleted from Helpers.jl with the typed input)
function get_dependent_value_with_ID(field_name::String,
                                     parameter::Dict,
                                     iID::Int64 = 1)
    return get_dependent_value(field_name, parameter)(iID)
end

function is_dependent(field_name::String, damage_parameter::Dict)
    if haskey(damage_parameter, field_name) && damage_parameter[field_name] isa Dict
        if !Data_Manager.has_key(damage_parameter[field_name]["Field"] * "NP1")
            @abort "$(damage_parameter[field_name]["Field"]) does not exist for value interpolation."
            return
        end
        field = Data_Manager.get_field(damage_parameter[field_name]["Field"], "NP1")
        return true, field
    end
    return false, nothing
end

function interpol_data(x::Union{Vector{Float64},Vector{Int64},Float64,Int64},
                       values::Dict{String,Any},
                       warning_flag::Bool = true)
    if warning_flag
        if values["min"] > minimum(x)
            @warn "Interpolation value is below interpolation range. Using minimum value of dataset."
        end
        if values["max"] < maximum(x)
            @warn "Interpolation value is above interpolation range. Using maximum value of dataset."
        end
        warning_flag = false
    end
    return evaluate(values["spl"], x)
end

abstract type AbstractDependentValue end

# Case 1: constant parameter (not field-dependent)
struct ConstantValue <: AbstractDependentValue
    value::Float64
end
(cv::ConstantValue)(iID::Int64) = cv.value

# Case 2: field-dependent, interpolated per node
struct InterpolatedValue{F,D} <: AbstractDependentValue
    field::F                      # e.g. Data_Manager field (NP1)
    data::D                       # parameter[field_name]["Data"]
    warning_flag::Base.RefValue{Bool}
end
function (iv::InterpolatedValue)(iID::Int64)
    val = interpol_data(iv.field[iID], iv.data, iv.warning_flag[])
    iv.warning_flag[] = false      # only warn once, on first call
    return val
end

"""
    get_dependent_value(field_name, parameter) -> AbstractDependentValue

Call once per field before a loop. Returns a callable `f(iID)` that
gives either the constant value or the interpolated field value.
"""
function get_dependent_value(field_name::String, parameter::Dict)
    dependent_value, dependent_field = is_dependent(field_name, parameter)
    if dependent_value
        return InterpolatedValue(dependent_field,
                                 parameter[field_name]["Data"],
                                 Ref(true))
    else
        return ConstantValue(parameter[field_name])
    end
end

function get_value(parameter::Union{Dict{Any,Any},Dict{String,Any}},
                   any_field_allocated::Bool,
                   key::String,
                   field_allocated::Bool)
    if field_allocated
        return Data_Manager.get_field(replace(key, " " => "_"))
    end
    if any_field_allocated
        if haskey(parameter, key)
            return Data_Manager.create_constant_node_scalar_field(replace(key, " " => "_"),
                                                                  Float64;
                                                                  default_value = parameter[key])
        else
            return Data_Manager.create_constant_node_scalar_field(replace(key, " " => "_"),
                                                                  Float64)
        end
    elseif haskey(parameter, key)
        return parameter[key]
    end

    return Float64(0.0)
end

function get_all_elastic_moduli(parameter::Union{Dict{Any,Any},Dict{String,Any}})
    state_factor_defined = haskey(parameter, "State Factor ID")

    if haskey(parameter, "Computed") &&
       !(state_factor_defined && Data_Manager.has_key("State Variables"))
        if parameter["Computed"]
            return
        end
    end

    bond_based = occursin("Bond-based", parameter["Material Model"])
    if bond_based
        bond_based = !occursin("Unified Bond-based", parameter["Material Model"])
    end
    bulk_field = Data_Manager.has_key("Bulk_Modulus")
    youngs_field = Data_Manager.has_key("Young's_Modulus")
    poissons_field = Data_Manager.has_key("Poisson's_Ratio")
    shear_field = Data_Manager.has_key("Shear_Modulus")

    any_field_allocated = bulk_field | youngs_field | poissons_field | shear_field |
                          state_factor_defined

    bulk = haskey(parameter, "Bulk Modulus") | bulk_field
    youngs = haskey(parameter, "Young's Modulus") | youngs_field
    shear = haskey(parameter, "Shear Modulus") | shear_field
    poissons = haskey(parameter, "Poisson's Ratio") | poissons_field

    K = get_value(parameter, any_field_allocated, "Bulk Modulus", bulk_field)
    E = get_value(parameter,
                  any_field_allocated,
                  "Young's Modulus",
                  youngs_field)
    G = get_value(parameter, any_field_allocated, "Shear Modulus", shear_field)

    nu = get_value(parameter,
                   any_field_allocated,
                   "Poisson's Ratio",
                   poissons_field)

    if bond_based
        nu_fixed = Data_Manager.get_dof() == 2 ? 1 / 3 : 1 / 4
        if nu != 0.0 && nu != nu_fixed
            @warn "Chosen Bond-based model only supports a fixed Poisson's ratio of " *
                  string(nu_fixed)
        end
        nu = nu_fixed
        poissons = true
    end
    if haskey(parameter, "Symmetry")
        symmetry = lowercase(parameter["Symmetry"])
        if occursin("anisotropic", symmetry)
            for iID in 1:6
                for jID in iID:6
                    if !haskey(parameter, "C" * string(iID) * string(jID))
                        @abort "C" * string(iID) * string(jID) * " not defined"
                        return
                    end
                end
            end
            return
        elseif occursin("transverse isotropic", symmetry)
            E_x = haskey(parameter, "Young's Modulus X")
            E_y = haskey(parameter, "Young's Modulus Y")
            nu_xy = haskey(parameter, "Poisson's Ratio XY")
            nu_yz = haskey(parameter, "Poisson's Ratio YZ")
            g_xy = haskey(parameter, "Shear Modulus XY")
            g_yz = haskey(parameter, "Shear Modulus YZ")
            if occursin("plane strain", symmetry)
                if !E_x || !E_y || !nu_xy || !nu_yz || !g_xy
                    @abort "Transverse isotropic material requires Young's Modulus X, Y, Poisson's Ratio XY, YZ, Shear Modulus XY"
                end
            elseif occursin("plane stress", symmetry)
                if !E_x || !E_y || !nu_xy || !g_xy
                    @abort "Transverse isotropic material requires Young's Modulus X, Y, Poisson's Ratio XY, Shear Modulus XY"
                end
            else
                if !E_x || !E_y || !nu_xy || !nu_yz || !g_xy || !g_yz
                    @abort "Transverse isotropic material requires Young's Modulus X, Y, Poisson's Ratio XY, YZ, Shear Modulus XY, YZ"
                end
            end
            return
        elseif occursin("orthotropic", symmetry)
            E_x = haskey(parameter, "Young's Modulus X")
            E_y = haskey(parameter, "Young's Modulus Y")
            E_z = haskey(parameter, "Young's Modulus Z")
            nu_xy = haskey(parameter, "Poisson's Ratio XY")
            nu_yz = haskey(parameter, "Poisson's Ratio YZ")
            nu_xz = haskey(parameter, "Poisson's Ratio XZ")
            g_xy = haskey(parameter, "Shear Modulus XY")
            g_yz = haskey(parameter, "Shear Modulus YZ")
            g_zx = haskey(parameter, "Shear Modulus XZ")
            if !E_x || !E_y || !E_z || !nu_xy || !nu_yz || !nu_xz || !g_xy || !g_yz || !g_zx
                @abort "Orthotropic material requires Young's Modulus X, Y, Z, Poisson's Ratio XY, YZ, XZ, Shear Modulus XY, YZ, XZ"
            end
            return
        end
    else
        @warn "Material symmetry is not defined, assuming isotropic material"
        parameter["Symmetry"] = "isotropic"
    end

    # tbd non isotropic material check
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
        state_factor = Data_Manager.get_field("State Variables")[:,
                                                                 parameter["State Factor ID"]]
        K .*= state_factor
        E .*= state_factor
        G .*= state_factor
    end

    parameter["Bulk Modulus"] = K
    parameter["Young's Modulus"] = E
    parameter["Shear Modulus"] = G
    parameter["Poisson's Ratio"] = nu
    parameter["Computed"] = true
    if any_field_allocated
        Data_Manager.get_field("Bulk_Modulus") .= K
        Data_Manager.get_field("Young's_Modulus") .= E
        Data_Manager.get_field("Shear_Modulus") .= G
        Data_Manager.get_field("Poisson's_Ratio") .= nu
    end
end

function get_symmetry(material::Dict)
    if !haskey(material, "Symmetry")
        return "3D"
    end
    if occursin("plane strain", lowercase(material["Symmetry"]))
        return "plane strain"
    end
    if occursin("plane stress", lowercase(material["Symmetry"]))
        return "plane stress"
    end
    return "3D"
end

function _dict_constant(parameter::Dict, key::String, id::Int64)
    if key in ("Poisson's Ratio", "Young's Modulus", "Shear Modulus")
        iID = parameter["Poisson's Ratio"] isa Float64 ? 1 : id
        return parameter[key][iID]
    end
    return get_dependent_value_with_ID(key, parameter, id)
end

function get_Hooke_matrix(parameter::Dict, symmetry::String, dof::Int64, ID::Int64 = 1)
    return _hooke_matrix((key, id) -> _dict_constant(parameter, key, id), symmetry, dof, ID)
end

function flaw_function(params::Dict,
                       coor::AbstractVector{<:Real},
                       stress::T)::Float64 where {T<:Union{Int64,Float64}}
    flaw = get(params, "Flaw Function", nothing)
    isnothing(flaw) && return stress

    if !haskey(flaw, "Active")
        @abort "Flaw Function needs an entry ''Active''."
    end
    if !haskey(flaw, "Function")
        @abort "Flaw Function needs an entry ''Function''."
    end

    flaw["Active"]::Bool || return stress

    if flaw["Function"] != "Pre-defined"
        @abort "Flaw Function ''$(flaw["Function"])'' is not implemented, " *
               "only ''Pre-defined'' is supported."
    end

    flaw_size::Float64 = flaw["Flaw Size"]
    flaw_magnitude::Float64 = flaw["Flaw Magnitude"]

    if !(0 < flaw_magnitude <= 1)
        @abort "Flaw Magnitude should be between 0 and 1"
    end
    if flaw_size <= 0
        @abort "Flaw Size must be positive."
    end

    # Squared distance without building a location vector: no allocation, no sqrt that
    # would only be squared again.
    dx = Float64(coor[1]) - Float64(flaw["Flaw Location X"])
    dy = Float64(coor[2]) - Float64(flaw["Flaw Location Y"])
    distance_squared = dx * dx + dy * dy

    if length(coor) == 3
        dz = Float64(coor[3]) - Float64(get(flaw, "Flaw Location Z", 0.0))
        distance_squared += dz * dz
    end

    return stress *
           (1 - flaw_magnitude * exp(-distance_squared / (flaw_size * flaw_size)))
end
end
