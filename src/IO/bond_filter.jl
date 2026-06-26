# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module Bond_Filter
using ...Parameter_Handling: get_bond_filters
using ....Data_Manager
using ....PeriLabExceptions: @abort
using DataFrames

include("../Core/Module_inclusion/set_Modules.jl")
global module_list = find_module_files(@__DIR__, "bond_filter_name")
for mod in module_list
    include(mod["File"])
end
export apply_bond_filters

"""
    apply_bond_filters(nlist::BondScalarState{Int64}, mesh::DataFrame, params::Dict, dof::Int64)

Apply the bond filters to the neighborhood list.

# Arguments
- `nlist::BondScalarState{Int64}`: The neighborhood list.
- `mesh::DataFrame`: The mesh.
- `params::Dict`: The parameters.
- `dof::Int64`: The degrees of freedom.
# Returns
- `nlist::BondScalarState{Int64}`: The filtered neighborhood list.
- `nlist_filtered_ids::BondScalarState{Int64}`: The filtered neighborhood list.
"""

function apply_bond_filters(nlist::BondScalarState{Int64},
                            mesh::DataFrame,
                            params::Dict,
                            dof::Int64)
    bond_filters = get_bond_filters(params)
    nlist_filtered_ids = nothing
    bond_norm = nothing
    contact_enabled = false
    if bond_filters[1]
        @debug "Apply bond filters"
        coor = names(mesh)[1:dof]
        nnodes = length(mesh[!, coor[1]])
        data = zeros(dof, nnodes)
        for i in 1:dof
            data[i, :] = values(mesh[!, coor[i]])
        end

        for (filter_name, filter) in bond_filters[2]
            contact_enabled = get(filter, "Allow Contact", false)
            if contact_enabled
                break
            end
        end
        if contact_enabled
            @info "Normal contact is applied within the bond filter."
            nlist_filtered_ids = fill(Vector{Int64}([]), nnodes)
            bond_norm = []
            for iID in 1:nnodes
                push!(bond_norm, [fill(1.0, dof) for n in 1:length(nlist[iID])])
            end
        end

        for (name, filter) in bond_filters[2]
            mod = create_module_specifics(filter["Type"],
                                          module_list,
                                          @__MODULE__,
                                          "bond_filter_name")
            if isnothing(mod)
                @warn "$(filter["Type"]) is not defined"
                return nlist, nlist_filtered_ids, bond_norm
            end
            filter_flag, normal = mod.run_bond_filter(nnodes, data, filter, nlist, dof)
            # Theoretically all bond filter can be in contact mode from memory side
            # but only the chosen ones are stored here.
            for iID in 1:nnodes
                if get(filter, "Allow Contact", false) &&
                   any(x -> x == false, filter_flag[iID])
                    indices = findall(x -> x in setdiff(nlist[iID],
                                                        nlist[iID][filter_flag[iID]]),
                                      nlist[iID])
                    nlist_filtered_ids[iID] = indices
                    for jID in indices
                        bond_norm[iID][jID] .= normal
                    end
                else
                    nlist[iID] = nlist[iID][filter_flag[iID]]
                end
            end
        end
        @debug "Finished applying bond filters"
    end
    return nlist, nlist_filtered_ids, bond_norm
end

end
