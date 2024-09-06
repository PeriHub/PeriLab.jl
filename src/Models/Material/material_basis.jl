# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using LinearAlgebra
using StaticArrays

function get_value(
    datamanager::Module,
    parameter::Union{Dict{Any,Any},Dict{String,Any}},
    any_field_allocated::Bool,
    key::String,
    field_allocated::Bool,
)
    if field_allocated
        return datamanager.get_field(replace(key, " " => "_"))
    end
    if any_field_allocated
        if haskey(parameter, key)
            return datamanager.create_constant_node_field(
                replace(key, " " => "_"),
                Float64,
                1,
                parameter[key],
            )
        else
            return datamanager.create_constant_node_field(
                replace(key, " " => "_"),
                Float64,
                1,
            )
        end
    elseif haskey(parameter, key)
        return parameter[key]
    end

    return Float64(0.0)
end

"""
    get_all_elastic_moduli(datamanager::Module, parameter::Union{Dict{Any,Any},Dict{String,Any}})

Returns the elastic moduli of the material.

# Arguments
- `parameter::Union{Dict{Any,Any},Dict{String,Any}}`: The material parameter.
"""
function get_all_elastic_moduli(
    datamanager::Module,
    parameter::Union{Dict{Any,Any},Dict{String,Any}},
)
    if haskey(parameter, "Computed")
        if parameter["Computed"]
            return nothing
        end
    end

    bulk_field = datamanager.has_key("Bulk_Modulus")
    Youngs_field = datamanager.has_key("Young's_Modulus")
    Poissons_field = datamanager.has_key("Poisson's_Ratio")
    shear_field = datamanager.has_key("Shear_Modulus")

    any_field_allocated = bulk_field | Youngs_field | Poissons_field | shear_field

    bulk = haskey(parameter, "Bulk Modulus") | bulk_field
    Youngs = haskey(parameter, "Young's Modulus") | Youngs_field
    Poissons = haskey(parameter, "Poisson's Ratio") | Poissons_field
    shear = haskey(parameter, "Shear Modulus") | shear_field

    K = get_value(datamanager, parameter, any_field_allocated, "Bulk Modulus", bulk_field)
    E = get_value(
        datamanager,
        parameter,
        any_field_allocated,
        "Young's Modulus",
        Youngs_field,
    )
    G = get_value(datamanager, parameter, any_field_allocated, "Shear Modulus", shear_field)
    nu = get_value(
        datamanager,
        parameter,
        any_field_allocated,
        "Poisson's Ratio",
        Poissons_field,
    )

    if bulk && Poissons
        E = 3 .* K .* (1 .- 2 .* nu)
        G = 3 .* K .* (1 .- 2 .* nu) ./ (2 .+ 2 .* nu)
    end
    if shear && Poissons
        E = 2 .* G .* (1 .+ nu)
        K = 2 .* G .* (1 .+ nu) ./ (3 .- 6 .* nu)
    end
    if bulk && shear
        E = 9 .* K .* G ./ (3 .* K .+ G)
        nu = (3 .* K .- 2 .* G) ./ (6 .* K .+ 2 .* G)
    end
    if Youngs && shear
        K = E .* G ./ (9 .* G .- 3 .* E)
        nu = E ./ (2 .* G) .- 1
    end

    if Youngs && bulk
        G = 3 .* K .* E ./ (9 .* K .- E)
        nu = (3 .* K .- E) ./ (6 .* K)
    end
    if Youngs && Poissons
        K = E ./ (3 .- 6 .* nu)
        G = E ./ (2 .+ 2 .* nu)
    end
    # tbd non isotropic material check
    if bulk + Youngs + shear + Poissons < 2
        @error "Minimum of two parameters are needed for isotropic material"
        return nothing
    end
    parameter["Bulk Modulus"] = K
    parameter["Young's Modulus"] = E
    parameter["Shear Modulus"] = G
    parameter["Poisson's Ratio"] = nu
    parameter["Computed"] = true
    if any_field_allocated
        datamanager.get_field("Bulk_Modulus") .= K
        datamanager.get_field("Young's_Modulus") .= E
        datamanager.get_field("Shear_Modulus") .= G
        datamanager.get_field("Poisson's_Ratio") .= nu
    end
end




"""
    get_Hooke_matrix(parameter::Dict, symmetry::String, dof::Int64, ID::Int64=1)

Returns the Hooke matrix of the material.

# Arguments
- `parameter::Union{Dict{Any,Any},Dict{String,Any}}`: The material parameter.
- `symmetry::String`: The symmetry of the material.
- `dof::Int64`: The degree of freedom.
- `ID::Int64=1`: ID of the point. Needed for point wise defined material properties.
# Returns
- `matrix::Matrix{Float64}`: The Hooke matrix.
"""
function get_Hooke_matrix(parameter::Dict, symmetry::String, dof::Int64, ID::Int64 = 1)
    """https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm"""

    if occursin("anisotropic", symmetry)
        aniso_matrix = @MMatrix zeros(Float64, 6, 6)
        for iID = 1:6
            for jID = iID:6
                if "C" * string(iID) * string(jID) in keys(parameter)
                    value = parameter["C"*string(iID)*string(jID)]
                else
                    @error "C" * string(iID) * string(jID) * " not defined"
                    return nothing
                end
                aniso_matrix[iID, jID] = value
                aniso_matrix[jID, iID] = value
            end
        end
        if dof == 3
            return aniso_matrix
        elseif occursin("plane strain", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            matrix[1:2, 1:2] = aniso_matrix[1:2, 1:2]
            matrix[3, 1:2] = aniso_matrix[6, 1:2]
            matrix[1:2, 3] = aniso_matrix[1:2, 6]
            matrix[3, 3] = aniso_matrix[6, 6]
            return matrix
        elseif occursin("plane stress", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            inv_aniso = inv(aniso_matrix)
            matrix[1:2, 1:2] = inv_aniso[1:2, 1:2]
            matrix[3, 1:2] = inv_aniso[6, 1:2]
            matrix[1:2, 3] = inv_aniso[1:2, 6]
            matrix[3, 3] = inv_aniso[6, 6]
            return inv(matrix)
        else
            @error "2D model defintion is missing; plane stress or plane strain "
            return nothing
        end
    end
    if !haskey(parameter, "Poisson's Ratio") ||
       !haskey(parameter, "Young's Modulus") ||
       !haskey(parameter, "Shear Modulus")
        @error "No valid definition of Hook matrix inputs."
        return nothing
    end
    iID = ID
    if parameter["Poisson's Ratio"] isa Float64
        iID = 1
    end
    if occursin("isotropic", symmetry)
        nu = parameter["Poisson's Ratio"][iID]
        E = parameter["Young's Modulus"][iID]
        G = parameter["Shear Modulus"][iID]
        temp = E / ((1 + nu) * (1 - 2 * nu))

        if dof == 3
            matrix = @MMatrix zeros(Float64, 2 * dof, 2 * dof)
            matrix[1, 1] = (1 - nu) * temp
            matrix[2, 2] = (1 - nu) * temp
            matrix[3, 3] = (1 - nu) * temp
            matrix[1, 2] = nu * temp
            matrix[2, 1] = nu * temp
            matrix[1, 3] = nu * temp
            matrix[3, 1] = nu * temp
            matrix[2, 3] = nu * temp
            matrix[3, 2] = nu * temp
            matrix[4, 4] = G
            matrix[5, 5] = G
            matrix[6, 6] = G
            return matrix
        elseif occursin("plane strain", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            matrix[1, 1] = (1 - nu) * temp
            matrix[2, 2] = (1 - nu) * temp
            matrix[3, 3] = G
            matrix[1, 2] = nu * temp
            matrix[2, 1] = nu * temp
            return matrix
        elseif occursin("plane stress", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            matrix[1, 1] = E / (1 - nu * nu)
            matrix[1, 2] = E * nu / (1 - nu * nu)
            matrix[2, 1] = E * nu / (1 - nu * nu)
            matrix[2, 2] = E / (1 - nu * nu)
            matrix[3, 3] = G
            return matrix
        else
            @error "2D model defintion is missing; plane stress or plane strain "
            return nothing
        end
    else
        matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)

        @warn "material model defintion is missing; assuming isotropic plane stress "
        nu = parameter["Poisson's Ratio"][iID]
        E = parameter["Young's Modulus"][iID]
        G = parameter["Shear Modulus"][iID]
        matrix[1, 1] = E / (1 - nu * nu)
        matrix[1, 2] = E * nu / (1 - nu * nu)
        matrix[2, 1] = E * nu / (1 - nu * nu)
        matrix[2, 2] = E / (1 - nu * nu)
        matrix[3, 3] = G
        return matrix

    end
end

"""
    distribute_forces(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, nlist_filtered_ids::SubArray, bond_force::SubArray, volume::SubArray, bond_damage::SubArray, displacements::SubArray, bond_norm::SubArray, force_densities::SubArray)

Distribute the forces on the nodes

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `nlist::SubArray`: The neighbor list.
- `nlist_filtered_ids::SubArray`:  The filtered neighbor list.
- `bond_force::SubArray`: The bond forces.
- `volume::SubArray`: The volumes.
- `bond_damage::SubArray`: The bond damage.
- `displacements::SubArray`: The displacements.
- `bond_norm::SubArray`: The pre defined bond normal.
- `force_densities::SubArray`: The force densities.
# Returns
- `force_densities::SubArray`: The force densities.
"""
function distribute_forces(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::SubArray,
    nlist_filtered_ids::SubArray,
    bond_force::SubArray,
    volume::SubArray,
    bond_damage::SubArray,
    displacements::SubArray,
    bond_norm::SubArray,
    force_densities::SubArray,
)

    for iID in nodes
        bond_mod = copy(bond_norm[iID])
        if length(nlist_filtered_ids[iID]) > 0
            for neighborID in nlist_filtered_ids[iID]
                if dot(
                    (displacements[nlist[iID][neighborID], :] - displacements[iID, :]),
                    bond_norm[iID][neighborID, :],
                ) > 0
                    bond_mod[neighborID, :] .= 0
                else
                    bond_mod[neighborID, :] = abs.(@view bond_norm[iID][neighborID, :])
                end
            end
        end

        force_densities[iID, :] .+= transpose(
            sum(
                bond_damage[iID] .* (@view bond_force[iID][:, :]) .* bond_mod .*
                volume[nlist[iID]],
                dims = 1,
            ),
        )

        # force_densities[nlist[iID], :] .-= bond_damage[iID] .* bond_force[iID] .* bond_mod .* volume[iID]
        force_densities[nlist[iID], :] .-=
            bond_damage[iID] .* (@view bond_force[iID][:, :]) .* (@view bond_mod[:, :]) .*
            volume[iID]
    end
    return force_densities
end

"""
    distribute_forces(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, bond_force::SubArray, volume::SubArray, bond_damage::SubArray, force_densities::SubArray)

Distribute the forces on the nodes

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `nlist::SubArray`: The neighbor list.
- `bond_force::SubArray`: The bond forces.
- `volume::SubArray`: The volumes.
- `bond_damage::SubArray`: The bond damage.
- `force_densities::SubArray`: The force densities.
# Returns
- `force_densities::SubArray`: The force densities.
"""
function distribute_forces(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::SubArray,
    bond_force::SubArray,
    volume::SubArray,
    bond_damage::SubArray,
    force_densities::SubArray,
)
    for iID in nodes
        force_densities[iID, :] .+= transpose(
            sum(bond_damage[iID] .* bond_force[iID] .* volume[nlist[iID]], dims = 1),
        )

        force_densities[nlist[iID], :] .-=
            bond_damage[iID] .* bond_force[iID] .* volume[iID]
    end
    return force_densities
end


"""
    matrix_to_voigt(matrix)

Convert a 2x2 or 3x3 matrix to Voigt notation (6x1 vector)

# Arguments
- `matrix::Matrix{Float64}`: The matrix.
# Returns
- `voigt::Vector{Float64}`: The Voigt notation.
"""
function matrix_to_voigt(matrix)
    if size(matrix) == (2, 2)
        return [matrix[1, 1]; matrix[2, 2]; 0.5 * (matrix[1, 2] + matrix[2, 1])]
    elseif size(matrix) == (3, 3)
        return [
            matrix[1, 1]
            matrix[2, 2]
            matrix[3, 3]
            0.5 * (matrix[2, 3] + matrix[3, 2])
            0.5 * (matrix[1, 3] + matrix[3, 1])
            0.5 * (matrix[1, 2] + matrix[2, 1])
        ]
    else
        @error "Unsupported matrix size for matrix_to_voigt"
        return nothing
    end
end

"""
    voigt_to_matrix(voigt::Union{Vector{Float64},Vector{Int64}})

Convert a Voigt notation (6x1 or 3x1 vector) to a 2x2 or 3x3 matrix

# Arguments
- `voigt::Vector{Float64}`: The Voigt notation.
# Returns
- `matrix::Matrix{Float64}`: The matrix.
"""
function voigt_to_matrix(voigt::Union{MVector,SVector,Vector})
    if length(voigt) == 3
        return [voigt[1] voigt[3]; voigt[3] voigt[2]]
    elseif length(voigt) == 6
        return [
            voigt[1] voigt[6] voigt[5]
            voigt[6] voigt[2] voigt[4]
            voigt[5] voigt[4] voigt[3]
        ]
    else
        @error "Unsupported matrix size for voigt_to_matrix"
        return nothing
    end
end

"""
    check_symmetry(prop::Dict, dof::Int64)

Check if the symmetry information is present in the material dictionary.

# Arguments
- `prop::Dict`: A dictionary containing material information.
- `dof::Int64`: The number of degrees of freedom.
# Returns
- `true`: If the symmetry information is present.
"""
function check_symmetry(prop::Dict, dof::Int64)
    if haskey(prop, "Symmetry")
        symmetry = prop["Symmetry"]
        if dof == 2
            if occursin("plane strain", symmetry) || occursin("plane stress", symmetry)
                return true
            else
                @error "Model definition is missing; plane stress or plane strain has to be defined for 2D"
                return
            end
        end
        return true
    end
end

"""
    flaw_function(params::Dict, coor::Union{Vector{Int64},Vector{Float64}}, stress::Float64)

Allows the modification of the yield stress at a specific position. This is typically used as starting point for plastic deformation.

# Arguments
- `params::Dict`: A dictionary containing material information.
- `coor::Union{Vector{Int64},Vector{Float64}, SubArray}`: Coordinate of the current point.
- `stress::Float64`: stresses to be modified.

# Returns
- `stress`::Float64: the modified stresses.
"""
function flaw_function(
    params::Dict,
    coor::Union{Vector{Int64},Vector{Float64},SubArray{Float64}},
    stress::Float64,
)
    if !haskey(params, "Flaw Function")
        return stress
    end
    if !haskey(params["Flaw Function"], "Active")
        @error "Flaw Function needs an entry ''Active''."
        return nothing
    end
    if !haskey(params["Flaw Function"], "Function")
        @error "Flaw Function needs an entry ''Function''."
        return nothing
    end
    if !params["Flaw Function"]["Active"]
        return stress
    end
    flaw_location::Vector{Float64} = zeros(length(coor))
    if params["Flaw Function"]["Function"] == "Pre-defined"
        flaw_size = params["Flaw Function"]["Flaw Size"]
        flaw_magnitude = params["Flaw Function"]["Flaw Magnitude"]
        if (1 < flaw_magnitude) || (flaw_magnitude < 0)
            @error "Flaw Magnitude should be between 0 and 1"
            return nothing
        end
        flaw_location[1] = params["Flaw Function"]["Flaw Location X"]
        flaw_location[2] = params["Flaw Function"]["Flaw Location Y"]
        if haskey(params["Flaw Function"], "Flaw location Z") && length(coor) == 3
            flaw_location[3] = params["Flaw Function"]["Flaw Location Z"]
        end
        modified_stress =
            stress * (
                1 -
                flaw_magnitude * exp(
                    -norm(coor - flaw_location) * norm(coor - flaw_location) / flaw_size /
                    flaw_size,
                )
            )

    else
        @warn "Not very user friendly right now"
        #global x = coor[1]
        #global y = coor[2]
        #
        #modified_stress = stress * (1 - eval(Meta.parse(params["Flaw Function"]["Function"])))
    end
    return modified_stress
end


"""
    get_symmetry(material::Dict)

Return the symmetry information from the given material dictionary.

# Arguments
- `material::Dict`: A dictionary containing material information.

# Returns
- If the key "Symmetry" is present in the dictionary, the corresponding value is returned.
- If the key is not present, the default value "3D" is returned.

# Example
```julia
material_dict = Dict("Symmetry" => "Cubic", "Color" => "Red")
symmetry = get_sym(material_dict)
"""
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

"""
    get_von_mises_stress(von_Mises_stress::Float64, dof::Int64, stress_NP1::Matrix{Float64})

# Arguments
- `von_Mises_stress::Float64`: Von Mises stress
- `dof::Int64`: Degree of freedom.
- `stress_NP1::Matrix{Float64}`: Stress.
# returns
- `spherical_stress_NP1::Float64`: Spherical stress
- `deviatoric_stress_NP1::Matrix{Float64}`: Deviatoric stress
"""
function get_von_mises_stress(
    von_Mises_stress::Float64,
    dof::Int64,
    stress_NP1::Matrix{Float64},
)

    spherical_stress_NP1 = sum(stress_NP1[i, i] for i = 1:dof) / 3
    deviatoric_stress_NP1 = stress_NP1[:, :] - spherical_stress_NP1 .* I(dof)

    von_Mises_stress =
        sqrt(3.0 / 2.0 * sum(deviatoric_stress_NP1[:, :] .* deviatoric_stress_NP1[:, :]))

    return von_Mises_stress, spherical_stress_NP1, deviatoric_stress_NP1
end

"""
    get_strain(stress_NP1::Matrix{Float64}, hooke_matrix::Matrix{Float64})

# Arguments
- `stress_NP1::Matrix{Float64}`: Stress.
- `hooke_matrix::Matrix{Float64}`: Hooke matrix
# returns
- `strain::Matrix{Float64}`: Strain
"""
function get_strain(stress_NP1::Matrix{Float64}, hooke_matrix::Matrix{Float64})
    return voigt_to_matrix(hooke_matrix' * matrix_to_voigt(stress_NP1))
end


function compute_Piola_Kirchhoff_stress(
    stress::Union{Matrix{Float64},SubArray{Float64}},
    deformation_gradient::Union{Matrix{Float64},SubArray{Float64}},
)
    return det(deformation_gradient) .* stress * invert(
        deformation_gradient,
        "Deformation gradient is singular and cannot be inverted.",
    )
end
