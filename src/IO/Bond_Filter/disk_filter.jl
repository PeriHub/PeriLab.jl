# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Disk_Filter
using LinearAlgebra
using .....Data_Manager
using .....PeriLabExceptions: @abort
export run_bond_filter, bond_filter_name
const TOLERANCE = 1.0e-14
"""
    bond_filter_name()

Return the name of this bond filter.

# Returns
- `String`: The name of the bond filter.
"""
function bond_filter_name()
    return "Disk"
end

"""
    bond_intersects_disc(p0::Vector{Float64}, p1::Vector{Float64}, center::Vector{Float64}, normal::Vector{Float64}, radius::Float64)

Check if a line segment intersects a disk.

# Arguments
- `p0::Vector{Float64}`: The start point of the line segment.
- `p1::Vector{Float64}`: The end point of the line segment.
- `center::Vector{Float64}`: The center of the disk.
- `normal::Vector{Float64}`: The normal of the plane.
- `radius::Float64`: The radius of the disk.
# Returns
- `Bool`: True if the line segment intersects the disk, False otherwise.
"""
function bond_intersects_disc(p0::Vector{Float64},
                              p1::Vector{Float64},
                              center::Vector{Float64},
                              normal::Vector{Float64},
                              radius::Float64)
    numerator = dot((center - p0), normal)
    denominator = dot((p1 - p0), normal)
    if abs(denominator) < TOLERANCE
        t = Inf
    else
        t = numerator / denominator
    end

    if t < 0.0 || t > 1.0
        return false
    end

    x = p0 + t .* (p1 - p0)
    distance = norm(x - center)

    if abs(distance) < radius^2
        return true
    end

    return false
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
    if dof == 3
        center = [filter["Center X"], filter["Center Y"], filter["Center Z"]]
        normal = [filter["Normal X"], filter["Normal Y"], filter["Normal Z"]]
    else
        @abort "Disk filter only implemented for 3D, use rectangular plane filter instead"
        return nothing
    end
    normal = normal ./ norm(normal)
    filter_flag::Vector{Vector{Bool}} = fill([], nnodes)
    for iID in 1:nnodes
        filter_flag[iID] = fill(true, length(nlist[iID]))
        for (jId, neighbor) in enumerate(nlist[iID])
            filter_flag[iID][jId] = !bond_intersects_disc(data[:, iID],
                                                          data[:, neighbor],
                                                          center,
                                                          normal,
                                                          filter["Radius"])
        end
    end
    return filter_flag, normal
end

end # module
