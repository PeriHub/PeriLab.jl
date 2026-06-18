# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Rectangular_Plane_Filter
using .....Data_Manager
export run_bond_filter, bond_filter_name

"""
    get_bond_filter_name()

Return the name of this bond filter.

# Returns
- `String`: The name of the bond filter.
"""
function bond_filter_name()
    return "Rectangular Plane"
end

"""
    bond_intersect_infinite_plane(p0::Vector{Float64}, p1::Vector{Float64}, lower_left_corner::Vector{Float64}, normal::Vector{Float64})

Check if a line segment intersects an infinite plane.

# Arguments
- `p0::Vector{Float64}`: The start point of the line segment.
- `p1::Vector{Float64}`: The end point of the line segment.
- `lower_left_corner::Vector{Float64}`: The lower left corner of the plane.
- `normal::Vector{Float64}`: The normal of the plane.
# Returns
- `Bool`: True if the line segment intersects the plane, False otherwise.
- `Vector{Float64}` or `undef`: The intersection point, or undef if no intersection.
"""
function bond_intersect_infinite_plane(p0::Vector{Float64},
                                       p1::Vector{Float64},
                                       lower_left_corner::Vector{Float64},
                                       normal::Vector{Float64})
    denominator = dot((p1 - p0), normal)
    if abs(denominator) < TOLERANCE
        return false, undef
    end
    t = dot((lower_left_corner - p0), normal) / denominator
    if 0.0 <= t <= 1.0
        return true, p0 + t .* (p1 - p0)
    end
    return false, undef
end

"""
    bond_intersect_rectangle_plane(x::Vector{Float64}, lower_left_corner::Vector{Float64}, bottom_unit_vector::Vector{Float64}, normal::Vector{Float64}, side_length::Float64, bottom_length::Float64)

Check if a point (intersection with the infinite plane) lies within the rectangle.

# Arguments
- `x::Vector{Float64}`: The point.
- `lower_left_corner::Vector{Float64}`: The lower left corner of the rectangle.
- `bottom_unit_vector::Vector{Float64}`: The unit vector along the bottom of the rectangle.
- `normal::Vector{Float64}`: The normal of the plane.
- `side_length::Float64`: The side length of the rectangle.
- `bottom_length::Float64`: The bottom length of the rectangle.
# Returns
- `Bool`: True if the point is inside the rectangle, False otherwise.
"""
function bond_intersect_rectangle_plane(x::Union{Vector{Float64},Vector{Int64}},
                                        lower_left_corner::Union{Vector{Float64},
                                                                 Vector{Int64}},
                                        bottom_unit_vector::Union{Vector{Float64},
                                                                  Vector{Int64}},
                                        normal::Union{Vector{Float64},Vector{Int64}},
                                        side_length::Real,
                                        bottom_length::Real)
    dr::Vector{Float64} = x - lower_left_corner
    bb::Float64 = dot(dr, bottom_unit_vector)
    if 0.0 <= bb && bb / bottom_length <= 1.0
        if length(normal) == 2
            return true
        end
        ua = cross(bottom_unit_vector, normal)
        aa = dot(dr, ua)
        if 0.0 <= aa && aa / side_length <= 1.0
            return true
        end
    end
    return false
end

"""
    run_bond_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::BondScalarState{Int64}, dof::Int64)

Apply the rectangular plane filter to the neighborhood list.

# Arguments
- `nnodes::Int64`: The number of nodes.
- `data::Matrix{Float64}`: The data.
- `filter::Dict`: The filter.
- `nlist::BondScalarState{Int64}`: The neighborhood list.
- `dof::Int64`: The degrees of freedom.
# Returns
- `filter_flag::Vector{Vector{Bool}}`: The filter flag.
- `normal::Vector{Float64}`: The normal vector of the plane.
"""
function run_bond_filter(nnodes::Int64,
                         data::Matrix{Float64},
                         filter::Dict,
                         nlist::BondScalarState{Int64},
                         dof::Int64)
    if dof == 2
        normal = [filter["Normal X"], filter["Normal Y"]]
        lower_left_corner = [filter["Lower Left Corner X"], filter["Lower Left Corner Y"]]
        bottom_unit_vector = [
            filter["Bottom Unit Vector X"],
            filter["Bottom Unit Vector Y"]
        ]
    else
        normal = [filter["Normal X"], filter["Normal Y"], filter["Normal Z"]]
        lower_left_corner = [
            filter["Lower Left Corner X"],
            filter["Lower Left Corner Y"],
            filter["Lower Left Corner Z"]
        ]
        bottom_unit_vector = [
            filter["Bottom Unit Vector X"],
            filter["Bottom Unit Vector Y"],
            filter["Bottom Unit Vector Z"]
        ]
    end
    normal = normal ./ norm(normal)
    bottom_unit_vector = bottom_unit_vector ./ norm(bottom_unit_vector)
    bottom_length = filter["Bottom Length"]
    side_length = filter["Side Length"]
    filter_flag::Vector{Vector{Bool}} = fill([], nnodes)
    for iID in 1:nnodes
        filter_flag[iID] = fill(true, length(nlist[iID]))
        for (jID, neighborID) in enumerate(nlist[iID])
            intersect_inf_plane,
            x = bond_intersect_infinite_plane(data[:, iID],
                                              data[:, neighborID],
                                              lower_left_corner,
                                              normal)
            bond_intersect = false
            if intersect_inf_plane
                bond_intersect = bond_intersect_rectangle_plane(x,
                                                                lower_left_corner,
                                                                bottom_unit_vector,
                                                                normal,
                                                                side_length,
                                                                bottom_length)
            end
            filter_flag[iID][jID] = !(intersect_inf_plane && bond_intersect)
        end
    end
    return filter_flag, normal
end

end # module
