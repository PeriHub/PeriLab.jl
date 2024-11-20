# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation
using DataStructures: OrderedDict
include("../../Support/geometry.jl")
using .Geometry
export pre_calculation_name
export init_model
export compute
export fields_for_local_synchronization

"""
    pre_calculation_name()

Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

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
    init_model(datamanager, nodes, parameter)

Inits the bond deformation calculation.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    parameter::Union{Dict,OrderedDict},
    block::Int64,
)
    dof = datamanager.get_dof()
    datamanager.create_bond_field("Deformed Bond Geometry", Float64, dof)
    datamanager.create_bond_field("Deformed Bond Length", Float64, 1)
    datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
    return datamanager
end

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}), parameter::Dict

Compute the bond deformation.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.
# Returns
- `datamanager`: Datamanager.
"""
function compute(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    parameter::Union{Dict,OrderedDict},
    block::Int64,
)
    nlist = datamanager.get_nlist()
    deformed_coor = datamanager.get_field("Deformed Coordinates", "NP1")
    deformed_coorN = datamanager.get_field("Deformed Coordinates", "N")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    Geometry.bond_geometry!(
        deformed_bond,
        deformed_bond_length,
        nodes,
        nlist,
        deformed_coor,
    )
    return datamanager
end



"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    # download_from_cores = false
    # upload_to_cores = true
    # datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end
end
