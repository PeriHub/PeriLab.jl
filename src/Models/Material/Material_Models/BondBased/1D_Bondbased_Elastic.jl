# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module OneD_Bond_based_elastic
using LoopVectorization
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
  init_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict{String, Any})

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    material_parameter::Dict{String,Any})
    constant = datamanager.create_constant_bond_field("Visual", Float64, 1)

    return datamanager
end

"""
    material_name()

Returns the name of the material model.
"""
function material_name()
    return "1D Bond-based Elastic"
end

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict{String, Any}, time::Float64, dt::Float64)

Calculate the elastic bond force for each node.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       material_parameter::Dict{String,Any},
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    undeformed_bond_length = datamanager.get_field("Bond Length")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_damage = datamanager.get_bond_damage("NP1")
    bond_force = datamanager.get_field("Bond Forces")
    coor = datamanager.get_field("Coordinates")
    E = material_parameter["Young's Modulus"]
    nlist = datamanager.get_nlist()
    # only one core
    id1 = material_parameter["Id1"]
    id2 = material_parameter["Id2"]
    idx = findfirst(==(id2), nlist[id1])
    bond_damage[id1][idx] = 0
    idx = findfirst(==(id1), nlist[id2])
    bond_damage[id2][idx] = 0
    c = 1.0
    for iID in nodes
        if any(deformed_bond_length[iID] .== 0)
            @error "Length of bond is zero due to its deformation."
            return nothing
        end

        # Calculate the bond force
        compute_bb_force!(bond_force[iID],
                          c,
                          bond_damage[iID],
                          deformed_bond_length[iID],
                          undeformed_bond_length[iID],
                          deformed_bond[iID])
        #bond_force[iID] =
        #    (
        #        0.5 .* constant[iID] .* bond_damage[iID] .*
        #        (deformed_bond_length[iID] .- undeformed_bond_length[iID]) ./
        #        undeformed_bond_length[iID]
        #    ) .* deformed_bond[iID] ./ deformed_bond_length[iID]
        #
    end
    # might be put in constant

    return datamanager
end

function compute_bb_force!(bond_force,
                           constant,
                           bond_damage,
                           deformed_bond_length,
                           undeformed_bond_length,
                           deformed_bond)
    for jID in eachindex(bond_force[:, 1])
        bond_force[jID, 1]
    end
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

end
