# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using ..Parameter_Handling: get_bond_filters

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
        # Line is parallel to the plane, may or may not lie on the plane
        # If it does lie on the plane, then the numerator will be zero
        # In either case, this function will return "no intersection"
        t = Inf
    else
        # The line intersects the plane
        t = numerator / denominator
    end

    if t < 0.0 || t > 1.0
        return false
    end

    # Intersection point
    x = p0 + t .* (p1 - p0)

    # Check if the intersection point is within the disk
    distance = norm(x - center)

    if abs(distance) < radius^2
        return true
    end

    return false
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
"""
function bond_intersect_infinite_plane(p0::Vector{Float64},
                                       p1::Vector{Float64},
                                       lower_left_corner::Vector{Float64},
                                       normal::Vector{Float64})
    denominator = dot((p1 - p0), normal)
    if abs(denominator) < TOLERANCE
        # Line is parallel to the plane
        # It may or may not lie on the plane
        # If it does lie on the plane, then the numerator will be zero
        # In either case, this function will return "no intersection"
        return false, undef
    end
    # The line intersects the plane

    t = dot((lower_left_corner - p0), normal) / denominator

    # Determine if the line segment intersects the plane
    if 0.0 <= t <= 1.0
        return true, p0 + t .* (p1 - p0)
    end
    # Intersection point
    return false, undef
end

"""
    bond_intersect_rectangle_plane(x::Vector{Float64}, lower_left_corner::Vector{Float64}, bottom_unit_vector::Vector{Float64}, normal::Vector{Float64}, side_length::Float64, bottom_length::Float64)

Check if a bond intersects a rectangle plane.

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
    if bond_filters[1]
        @debug "Apply bond filters"
        coor = names(mesh)[1:dof]
        nnodes = length(mesh[!, coor[1]])
        data = zeros(dof, nnodes)
        for i in 1:dof
            data[i, :] = values(mesh[!, coor[i]])
        end
        #TODO to the bottom, because right now all filters have contact if true
        contact_enabled = false
        for (name, filter) in bond_filters[2]
            contact_enabled = get(filter, "Allow Contact", false)
            if contact_enabled
                break
            end
        end
        if contact_enabled
            nlist_filtered_ids = fill(Vector{Int64}([]), nnodes)
            bond_norm = []
            for iID in 1:nnodes
                push!(bond_norm, [fill(1.0, dof) for n in 1:length(nlist[iID])])
            end
        end

        for (name, filter) in bond_filters[2]
            if filter["Type"] == "Disk"
                filter_flag, normal = disk_filter(nnodes, data, filter, nlist, dof)
            elseif filter["Type"] == "Rectangular_Plane"
                filter_flag,
                normal = rectangular_plane_filter(nnodes, data, filter, nlist,
                                                  dof)
            end
            for iID in 1:nnodes
                if contact_enabled && any(x -> x == false, filter_flag[iID])
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

"""
    disk_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::BondScalarState{Int64}, dof::Int64)

Apply the disk filter to the neighborhood list.

# Arguments
- `nnodes::Int64`: The number of nodes.
- `data::Matrix{Float64}`: The data.
- `filter::Dict`: The filter.
- `nlist::BondScalarState{Int64}`: The neighborhood list.
- `dof::Int64`: The degrees of freedom.
# Returns
- `filter_flag::Vector{Bool}`: The filter flag.
- `normal::Vector{Float64}`: The normal vector of the disk.
"""
function disk_filter(nnodes::Int64,
                     data::Matrix{Float64},
                     filter::Dict,
                     nlist::BondScalarState{Int64},
                     dof::Int64)
    if dof == 3
        center = [filter["Center X"], filter["Center Y"], filter["Center Z"]]
        normal = [filter["Normal X"], filter["Normal Y"], filter["Normal Z"]]
    else
        @error "Disk filter only implemented for 3D, use rectangular plane filter instead"
        return nothing
    end
    #normalize vector
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

"""
    rectangular_plane_filter(nnodes::Int64, data::Matrix{Float64}, filter::Dict, nlist::BondScalarState{Int64}, dof::Int64)

Apply the rectangular plane filter to the neighborhood list.

# Arguments
- `nnodes::Int64`: The number of nodes.
- `data::Matrix{Float64}`: The data.
- `filter::Dict`: The filter.
- `nlist::BondScalarState{Int64}`: The neighborhood list.
- `dof::Int64`: The degrees of freedom.
# Returns
- `filter_flag::Vector{Bool}`: The filter flag.
- `normal::Vector{Float64}`: The normal vector of the disk.
"""
function rectangular_plane_filter(nnodes::Int64,
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
    #normalize vector
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
