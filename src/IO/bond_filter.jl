# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module Bond_Filter
using ....InputDeck: BondFilterParams, bond_filter
using ....Data_Manager
using ....ModuleLoader: find_registered_modules
using DataFrames
using LinearAlgebra: dot

for file in find_registered_modules(joinpath(@__DIR__, "Bond_Filter"),
                                    "register_bond_filter")
    include(file)
end
export apply_bond_filters

"""
    apply_bond_filters(nlist::NestedNeighborhoodList, mesh::DataFrame, filters::Dict{String,BondFilterParams}, dof::Int64)

Apply the bond filters to the neighborhood list. Filters without `Allow Contact` remove
the bonds; with it the bonds stay as contact bonds (see `Material.distribute_force_densities`).

# Arguments
- `nlist::NestedNeighborhoodList`: The neighborhood list.
- `mesh::DataFrame`: The mesh.
- `filters::Dict{String,BondFilterParams}`: The bond filters of the input deck.
- `dof::Int64`: The degrees of freedom.
# Returns
- `nlist::NestedNeighborhoodList`: The filtered neighborhood list.
- `bond_norm`: `nothing` without contact filters; else a vector per bond: zero for a
  normal bond, for a contact bond the filter normal oriented along the bond and scaled by
  the `Contact Stiffness Factor`.
"""
function apply_bond_filters(nlist::NestedNeighborhoodList,
                            mesh::DataFrame,
                            filters::Dict{String,BondFilterParams},
                            dof::Int64)
    bond_norm = nothing
    isempty(filters) && return nlist, bond_norm
    @debug "Apply bond filters"
    data = permutedims(Matrix{Float64}(mesh[!, names(mesh)[1:dof]]))
    nnodes = size(data, 2)
    if any(filter -> filter.allow_contact, values(filters))
        @info "Normal contact is applied within the bond filter."
        bond_norm = [[zeros(dof) for _ in neighbors] for neighbors in nlist]
    end
    for (name, filter) in filters
        # the type was checked when the input deck was read
        filter_flag, normal = bond_filter(filter.type).run(nnodes, data, filter, nlist,
                                                           dof)
        for iID in 1:nnodes
            if filter.allow_contact
                for jID in findall(!, filter_flag[iID])
                    bond = data[:, nlist[iID][jID]] - data[:, iID]
                    bond_norm[iID][jID] .= sign(dot(bond, normal)) *
                                           filter.contact_stiffness_factor .* normal
                end
            else
                nlist[iID] = nlist[iID][filter_flag[iID]]
                isnothing(bond_norm) || (bond_norm[iID] = bond_norm[iID][filter_flag[iID]])
            end
        end
    end
    @debug "Finished applying bond filters"
    return nlist, bond_norm
end

end
