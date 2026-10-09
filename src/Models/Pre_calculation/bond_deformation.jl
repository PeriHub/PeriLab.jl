# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation
using .......Data_Manager
using .......Geometry: bond_geometry!
export pre_calculation_name
export init_model
export compute
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_pre_calculation
"Switch only: this pre-calculation has no parameters."
@params struct DeformedBondGeometryParams
end
__init__() = register_pre_calculation("Deformed Bond Geometry", DeformedBondGeometryParams)

"""
    pre_calculation_name()

Gives the pre_calculation name. PeriLab loads the module because it defines this function; the input deck uses the name passed to `register_*` in `__init__()`.

# Arguments

# Returns
- `name::String`: The name of the Pre_Calculation.

Example:
```julia
println(pre_calculation_name())
"Deformed Bond Geometry"
```
"""
function pre_calculation_name()
    return "Deformed Bond Geometry"
end

"""
    init_model(nodes, block)

Inits the bond deformation calculation.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.

"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    dof = Data_Manager.get_dof()
    Data_Manager.create_bond_vector_state("Deformed Bond Geometry", Float64, dof)
    Data_Manager.create_bond_scalar_state("Deformed Bond Length", Float64)
    Data_Manager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
end

"""
    compute(nodes::AbstractVector{Int64}, block::Int64)

Compute the bond deformation.

# Arguments
- `nodes`: List of nodes.
"""
function compute(nodes::AbstractVector{Int64}, block::Int64)
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    deformed_coor::NodeVectorField{Float64} = Data_Manager.get_field("Deformed Coordinates",
                                                                     "NP1")
    deformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Deformed Bond Geometry",
                                                                     "NP1")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")
    bond_geometry!(deformed_bond,
                   deformed_bond_length,
                   nodes,
                   nlist,
                   deformed_coor)
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    # download_from_cores = false
    # upload_to_cores = true
    # Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return
end
end
