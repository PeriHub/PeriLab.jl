# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Flow
using LinearAlgebra
using StaticArrays
using TimerOutputs: @timeit
using .....Data_Manager
using .....PeriLabExceptions: @abort
using .....Helpers: node_tensor
export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_thermal

@enum ThermalFlowType BondBased Correspondence

@params "Heat conduction between nodes" struct ThermalFlowParams
    "Formulation of the heat flow: per bond, or from the temperature gradient as in correspondence materials"
    type::ThermalFlowType = opt("Type"; default = BondBased)
end

__init__() = register_thermal("Thermal Flow", ThermalFlowParams)

"""
    init_model(nodes, p, thermal, block)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::ThermalFlowParams`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block (`thermal.base`: Thermal Conductivity).
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::ThermalFlowParams, thermal,
                    block::Int64)
    if thermal.base.thermal_conductivity === nothing
        @abort "Thermal Conductivity not defined."
    end
    p.type == BondBased || return nothing
    weight = Data_Manager.create_constant_bond_scalar_state("Thermal Bond Weight", Float64)
    init_bond_weights!(weight, nodes, Data_Manager.get_dof(), Data_Manager.get_nlist(),
                       first(thermal.base.thermal_conductivity),
                       Data_Manager.get_field("Bond Length"),
                       Data_Manager.get_field("Horizon"), Data_Manager.get_field("Volume"))
end

"""
    init_bond_weights!(weight, nodes, dof, nlist, lambda, undeformed_bond_length, horizon, volume)

Constant part of the bond-based heat flow: `lambda kernel V_j / |xi|` with the kernel
`6 / (pi delta^(dof+1))` [OterkusS2014b](@cite).
"""
function init_bond_weights!(weight, nodes::AbstractVector{Int64}, dof::Int64, nlist,
                            lambda::Float64, undeformed_bond_length, horizon, volume)
    for iID in nodes
        kernel = 6.0 / (pi * horizon[iID]^(dof + 1))
        for (jID, neighborID) in enumerate(nlist[iID])
            weight[iID][jID] = lambda * kernel * volume[neighborID] /
                               undeformed_bond_length[iID][jID]
        end
    end
end

"""
    compute_model(nodes, p, thermal, block, time, dt)

Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::ThermalFlowParams`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block (`thermal.base`: Thermal Conductivity).
- `block::Int64`: Current block
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64}, p::ThermalFlowParams, thermal,
                       block::Int64, time::Float64, dt::Float64)
    dof::Int64 = Data_Manager.get_dof()
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    heat_flow::NodeScalarField{Float64} = Data_Manager.get_field("Heat Flow", "NP1")

    undeformed_bond::BondVectorField{Float64} = Data_Manager.get_field("Bond Geometry")
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
    temperature::NodeScalarField{Float64} = Data_Manager.get_field("Temperature", "NP1")
    active = Data_Manager.get_field("Active")

    lambda::Float64 = thermal.base.thermal_conductivity
    rotation::Bool = Data_Manager.get_element_rotation()
    rotation_tensor = nothing
    if rotation
        rotation_tensor = Data_Manager.get_field("Rotation Tensor")::Array{Float64,3}
    end
    if p.type == BondBased
        @timeit "heat_flow" compute_heat_flow_state_bond_based!(nodes, nlist,
                                                                Data_Manager.get_field("Thermal Bond Weight"),
                                                                bond_damage, active,
                                                                temperature, heat_flow)
        return

    elseif p.type == Correspondence
        lambda_matrix = @MMatrix zeros(Float64, dof, dof)
        Kinv = Data_Manager.get_field("Inverse Shape Tensor")
        if length(lambda) == 1
            for i in 1:dof
                lambda_matrix[i, i] = lambda
            end
        else
            for i in 1:dof
                lambda_matrix[i, i] = lambda[i]
            end
        end
        @timeit "heat_flow_correspondence" compute_heat_flow_state_correspondence!(nodes,
                                                                                   dof,
                                                                                   nlist,
                                                                                   lambda_matrix,
                                                                                   rotation_tensor,
                                                                                   bond_damage,
                                                                                   undeformed_bond,
                                                                                   Kinv,
                                                                                   temperature,
                                                                                   volume,
                                                                                   heat_flow)
    end
end

"""
[BrighentiR2021](@cite)
is a prototype with some errors
"""
function compute_heat_flow_state_correspondence!(nodes::AbstractVector{Int64},
                                                 dof::Int64,
                                                 nlist::BondScalarState{Int64},
                                                 lambda::Union{Matrix{Float64},MMatrix},
                                                 rotation_tensor::Union{Nothing,
                                                                        Array{Float64,3}},
                                                 bond_damage::BondScalarState{Float64},
                                                 undeformed_bond::BondVectorState{Float64,
                                                                                  D},
                                                 Kinv::Array{Float64,3},
                                                 temperature::NodeScalarField{Float64},
                                                 volume::NodeScalarField{Float64},
                                                 heat_flow::NodeScalarField{Float64}) where {D}
    conductivity = SMatrix{D,D,Float64,D * D}(lambda)
    @inbounds for iID in nodes
        neighbors = nlist[iID]
        bonds = undeformed_bond[iID]
        damage = bond_damage[iID]
        H = zero(SVector{D,Float64})
        for jID in eachindex(neighbors)
            neighborID = neighbors[jID]
            temp_state = (temperature[neighborID] - temperature[iID]) *
                         volume[neighborID] *
                         damage[jID]
            H += temp_state * bonds[jID]
        end
        K = node_tensor(Kinv, iID, Val(D))
        nablaT = K * H
        q = rotated_conductivity(conductivity, rotation_tensor, iID) * nablaT
        for jID in eachindex(neighbors)
            neighborID = neighbors[jID]
            flow = dot(K * bonds[jID], q)
            heat_flow[iID] -= flow * volume[neighborID]
            heat_flow[neighborID] += flow * volume[iID]
        end
    end
end

# the conductivity in the global frame: R Λ Rᵀ (as rotate_second_order_tensor2x2/3x3)
rotated_conductivity(conductivity, ::Nothing, iID::Int64) = conductivity
function rotated_conductivity(conductivity::SMatrix{D,D,Float64},
                              rotation_tensor::Array{Float64,3}, iID::Int64) where {D}
    R = node_tensor(rotation_tensor, iID, Val(D))
    return R * conductivity * R'
end

"""
    compute_heat_flow_state_bond_based!(nodes, nlist, weight, bond_damage, active, temperature, heat_flow)

Bond-based heat flow [OterkusS2014b](@cite): `heat_flow_i -= sum_j weight_ij d_ij (T_j - T_i)`
with the bond damage `d_ij` and the constant bond weight of `init_bond_weights!`.
"""
function compute_heat_flow_state_bond_based!(nodes::AbstractVector{Int64},
                                             nlist::BondScalarState{Int64},
                                             weight::BondScalarState{Float64},
                                             bond_damage::BondScalarState{Float64},
                                             active::AbstractVector{Bool},
                                             temperature::AbstractVector{Float64},
                                             heat_flow::AbstractVector{Float64})
    @inbounds for iID in nodes
        active[iID] || continue
        Ti = temperature[iID]
        flow = 0.0
        w = weight[iID]
        d = bond_damage[iID]
        for (jID, neighborID) in enumerate(nlist[iID])
            flow += w[jID] * d[jID] * (temperature[neighborID] - Ti)
        end
        heat_flow[iID] -= flow
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
