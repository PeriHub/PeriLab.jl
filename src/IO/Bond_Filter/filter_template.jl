# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Filter_template
using .....Data_Manager
export run_bond_filter, bond_filter_name
const TOLERANCE = 1.0e-14
"""
    bond_filter_name()

Return the name of this bond filter.

# Returns
- `String`: The name of the bond filter.
"""
function bond_filter_name()
    return "Template"
end

"""
    run_bond_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::BondScalarState{Int64}, dof::Int64)

Apply the disk filter to the neighborhood list.

# Arguments
- `nnodes::Int64`: The number of nodes.
- `data::Matrix{Float64}`: The data.
- `filter::Dict`: The filter.
- `nlist::BondScalarState{Int64}`: The neighborhood list.
- `dof::Int64`: The degrees of freedom.
# Returns
- `filter_flag::Vector{Vector{Bool}}`: The filter flag.
- `normal::Vector{Float64}`: The normal vector of the disk.
"""
function run_bond_filter(nnodes::Int64,
                         data::Matrix{Float64},
                         filter::Dict,
                         nlist::BondScalarState{Int64},
                         dof::Int64)
    @info "please add your filter here"
end

end
