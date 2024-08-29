# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Deformation_Gradient
include("../../Support/geometry.jl")
using .Geometry: compute_deformation_gradient
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
"Deformation Gradient"
```
"""
function pre_calculation_name()
    return "Deformation Gradient"
end

"""
    init_model(datamanager, nodes, corrosion_parameter)

Inits the bond-based corrosion model. This template has to be copied, the file renamed and edited by the user to create a new corrosion. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `corrosion parameter::Dict(String, Any)`: Dictionary with corrosion parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    corrosion_parameter::Dict,
    block::Int64,
)
    dof = datamanager.get_dof()
    datamanager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    return datamanager
end

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Compute the deformation gradient.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
# Returns
- `datamanager`: Datamanager.
"""
function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")

    deformation_gradient = compute_deformation_gradient(
        nodes,
        dof,
        nlist,
        volume,
        omega,
        bond_damage,
        deformed_bond,
        undeformed_bond,
        inverse_shape_tensor,
        deformation_gradient,
    )

    return datamanager
end

end
