# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Shape_Tensor
using DataStructures: OrderedDict
include("../../Support/geometry.jl")
using .Geometry
export pre_calculation_name
export init_model
export compute


"""
    pre_calculation_name()

Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the Pre_Calculation.

Example:
```julia
println(pre_calculation_name())
"Shape Tensor"
```
"""
function pre_calculation_name()
    return "Shape Tensor"
end


"""
    init_model(datamanager, nodes, parameter)

Inits the shape tensor calculation.

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
)
    dof = datamanager.get_dof()
    datamanager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
    datamanager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)
    return datamanager
end

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Compute the shape tensor.

# Arguments
- `datamanager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
# Returns
- `datamanager`: Datamanager.
"""
function compute(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    parameter::Union{Dict,OrderedDict},
)
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    shape_tensor = datamanager.get_field("Shape Tensor")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")
    shape_tensor, inverse_shape_tensor = Geometry.shape_tensor(
        nodes,
        dof,
        nlist,
        volume,
        omega,
        bond_damage,
        undeformed_bond,
        shape_tensor,
        inverse_shape_tensor,
    )
    return datamanager
end



end
