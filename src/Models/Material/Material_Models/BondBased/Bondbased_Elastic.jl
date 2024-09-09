# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Elastic
include("../../material_basis.jl")
using TimerOutputs
export init_model
export fe_support
export material_name
export compute_model

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
  init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    material_parameter::Dict,
)
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
    fields_for_local_synchronization()

Returns a user developer defined local synchronization. This happens before each model.

The structure of the Dict must because

    synchfield = Dict(
        "Field name" =>
            Dict("upload_to_cores" => true, "dof" => datamanager.get_dof()),
    )

or

    synchfield = Dict(
        "Field name" =>
            Dict("download_from_cores" => true, "dof" => datamanager.get_dof()),
    )

# Arguments

"""
function fields_for_local_synchronization()
    return Dict()
end

"""
    compute_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)

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
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    material_parameter::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)
    # global dof
    # global horizon
    dof = datamanager.get_dof()
    horizon = datamanager.get_field("Horizon")
    symmetry::String = get_symmetry(material_parameter)
    undeformed_bond_length = datamanager.get_field("Bond Length")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_damage = datamanager.get_bond_damage("NP1")
    bond_force = datamanager.get_field("Bond Forces")

    E = material_parameter["Young's Modulus"]

    for iID in nodes
        if symmetry == "plane stress"
            constant = 12 / (2 * (1 - 1 / 3)) / (pi * horizon[iID]^3) # from EQ 2.9 +2.9 D=2 in Handbook of PD
        elseif symmetry == "plane strain"
            constant = 12 / (2 * (1 - 0.25 + 0.25 * 0.25)) / (pi * horizon[iID]^3) # from EQ 2.12 + 2.9 D=2 in Handbook of PD
        else
            constant = 18 / (3 - 2 / 3) / (pi * horizon[iID]^4) # from EQ 2.12 D=3 in Handbook of PD
        end
        if any(deformed_bond_length[iID] .== 0)
            @error "Length of bond is zero due to its deformation."
            return nothing
        end
        # Calculate the bond force
        bond_force[iID] =
            (
                0.5 .* constant .* bond_damage[iID] .*
                (deformed_bond_length[iID] .- undeformed_bond_length[iID]) ./
                undeformed_bond_length[iID]
            ) .* deformed_bond[iID] ./ deformed_bond_length[iID]

    end
    # checks if E is scalar or a vector. Is needed for point wise definition
    #bond_force[nodes] .*= isa(E, Float64) ? E : E[nodes]
    if isa(E, Float64) # faster than the on line solution
        bond_force[nodes] .*= E
    else
        bond_force[nodes] .*= E[nodes]
    end
    return datamanager
end

end
