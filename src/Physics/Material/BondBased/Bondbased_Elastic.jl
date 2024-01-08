# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Elastic
include("../material_basis.jl")
export init_material_model
export fe_support
export init_material_model
export material_name
export compute_forces

"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
    return false
end

"""
  init_material_model(datamanager::Module)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module)
    return datamanager
end


"""
    material_name()

Returns the name of the material model.
"""
function material_name()
    return "Bond-based Elastic"
end

"""
    compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)

Calculate the elastic bond force for each node.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64, to::TimerOutput)
    # global dof
    # global horizon
    dof = datamanager.get_dof()
    horizon = datamanager.get_field("Horizon")
    symmetry::String = get_symmmetry(material_parameter)
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    bond_force = datamanager.get_field("Bond Forces")

    K = material_parameter["Bulk Modulus"]
    E = material_parameter["Young's Modulus"]
    for iID in nodes
        if symmetry == "plane stress"
            constant = 12.0 * E / (pi * (1 + 1.0 / 3) * horizon[iID]^3)
        elseif symmetry == "plane strain"
            constant = 12.0 * E / (pi * (1.25) * horizon[iID]^3)
        else
            constant = 18.0 * K / (pi * horizon[iID]^4)
        end
        if any(deformed_bond[iID][:, end] .== 0)
            @error "Length of bond is zero due to its deformation."
            return nothing
        end
        # Calculate the bond force
        bond_force[iID] = (0.5 .* constant .* bond_damage[iID][:] .* (deformed_bond[iID][:, end] .- undeformed_bond[iID][:, end]) ./ undeformed_bond[iID][:, end]) .* deformed_bond[iID][:, 1:dof] ./ deformed_bond[iID][:, end]

    end

    return datamanager
end

end