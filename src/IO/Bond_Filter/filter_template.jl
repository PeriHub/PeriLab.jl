# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Filter_template
using .....Data_Manager
using .....InputDeck: BondFilterParams, register_bond_filter
export run_bond_filter
const TOLERANCE = 1.0e-14

# Uncomment to make the filter available as `Type: Template`; `required` lists
# the Bond Filters keys it needs.
# __init__() = register_bond_filter("Template", run_bond_filter; required = String[])

"""
    run_bond_filter(nnodes::Int64, data::Matrix{Float64}, filter::BondFilterParams, nlist::BondScalarState{Int64}, dof::Int64)

Apply the disk filter to the neighborhood list.

# Arguments
- `nnodes::Int64`: The number of nodes.
- `data::Matrix{Float64}`: The data.
- `filter::BondFilterParams`: The filter.
- `nlist::BondScalarState{Int64}`: The neighborhood list.
- `dof::Int64`: The degrees of freedom.
# Returns
- `filter_flag::Vector{Vector{Bool}}`: The filter flag.
- `normal::Vector{Float64}`: The normal vector of the disk.
"""
function run_bond_filter(nnodes::Int64,
                         data::Matrix{Float64},
                         filter::BondFilterParams,
                         nlist::BondScalarState{Int64},
                         dof::Int64)
    @info "please add your filter here"
end

end
