# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module OneD_Bond_Based_Elastic
using LoopVectorization

using .......Data_Manager
using ......ParameterSpec: @params, register_material
using .......PeriLabExceptions: @abort

export init_model
export fe_support
export compute_model

@params struct OneDBondbasedElasticParams
    id1::Int64 = req("Id1"; min = 1)
    id2::Int64 = req("Id2"; min = 1)
end
__init__() = register_material("1D Bond-based Elastic", OneDBondbasedElasticParams)

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
  init_model(nodes::AbstractVector{Int64}, p, material, block::Int64)

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `p`: Model parameters; `material::BlockMaterial`: typed block material (`material.base`, `material.moduli`, `material.symmetry`).
"""
function init_model(nodes::AbstractVector{Int64},
                    p::OneDBondbasedElasticParams,
                    material,
                    block::Int64)
    constant = Data_Manager.create_constant_bond_scalar_state("Visual", Float64)
end

"""
    compute_model(nodes::AbstractVector{Int64}, p, material, time::Float64, dt::Float64)

Calculate the elastic bond force for each node.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p`: Model parameters; `material::BlockMaterial`: typed block material (`material.base`, `material.moduli`, `material.symmetry`).
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64},
                       p::OneDBondbasedElasticParams,
                       material,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    undeformed_bond_length = Data_Manager.get_field("Bond Length")
    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    deformed_bond = Data_Manager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = Data_Manager.get_field("Deformed Bond Length", "NP1")
    bond_damage = Data_Manager.get_bond_damage("NP1")
    bond_force = Data_Manager.get_field("Bond Forces")
    coor = Data_Manager.get_field("Coordinates")
    E = material.moduli.youngs_modulus
    nlist = Data_Manager.get_nlist()
    # only one core
    id1 = p.id1
    id2 = p.id2
    idx = findfirst(==(id2), nlist[id1])
    bond_damage[id1][idx] = 0
    idx = findfirst(==(id1), nlist[id2])
    bond_damage[id2][idx] = 0
    c = 1.0
    for iID in nodes
        if any(deformed_bond_length[iID] .== 0)
            @abort "Length of bond is zero due to its deformation."
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
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
