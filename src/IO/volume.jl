# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    calculate_volume(element_type::String, vertices::Vector{Vector{Float64}})

Calculate the volume of a element.

# Arguments
- `element_type`: The element type of the element.
- `vertices`: The vertices of the element.
# Returns
- `volume`: The volume of the element.
"""
function calculate_volume(element_type::String, vertices::Vector{Vector{Float64}})
    if element_type in ["Quad4", "CPS3"]
        return area_of_polygon(vertices)
    elseif element_type == "Tet4"
        return tetrahedron_volume(vertices)
    elseif element_type == "Wedge6"
        return wedge6_volume(vertices)
    elseif element_type == "Hex8"
        return hex8_volume(vertices)
    else
        @error "Element type $element_type currently not supported"
        return nothing
    end
end

"""
    tetrahedron_volume(tet_vertices)

Calculate the volume of a tetrahedron.

# Arguments
- `tet_vertices`: The vertices of the tetrahedron.
# Returns
- `volume`: The volume of the tetrahedron.
"""
function tetrahedron_volume(tet_vertices::Vector{Vector{Float64}})
    mat = hcat(hcat(tet_vertices...)', ones(4))  # Augmenting matrix with ones in the fourth column
    volume = abs(det(mat) / 6)   # Using det function to calculate determinant
    return volume
end

"""
hex8_volume(hex_vertices)

Calculate the volume of a hex.

# Arguments
- `hex_vertices`: The vertices of the wedge.
# Returns
- `volume`: The volume of the wedge.
"""
function hex8_volume(hex_vertices::Vector{Vector{Float64}})
    tets = [
        [hex_vertices[1], hex_vertices[2], hex_vertices[4], hex_vertices[5]],
        [hex_vertices[2], hex_vertices[3], hex_vertices[4], hex_vertices[7]],
        [hex_vertices[2], hex_vertices[5], hex_vertices[6], hex_vertices[7]],
        [hex_vertices[4], hex_vertices[5], hex_vertices[7], hex_vertices[8]],
        [hex_vertices[2], hex_vertices[4], hex_vertices[5], hex_vertices[7]]
    ]

    volume = 0
    for tet in tets
        volume += tetrahedron_volume(tet)
    end

    return volume
end

"""
hex8_volume(hex_vertices)

Calculate the volume of a hex.

# Arguments
- `hex_vertices`: The vertices of the wedge.
# Returns
- `volume`: The volume of the wedge.
"""
function hex8_volume(hex_vertices::Matrix{Float64})
    tets = [
        [hex_vertices[:, 1], hex_vertices[:, 2], hex_vertices[:, 4], hex_vertices[:, 5]],
        [hex_vertices[:, 2], hex_vertices[:, 3], hex_vertices[:, 4], hex_vertices[:, 7]],
        [hex_vertices[:, 2], hex_vertices[:, 5], hex_vertices[:, 6], hex_vertices[:, 7]],
        [hex_vertices[:, 4], hex_vertices[:, 5], hex_vertices[:, 7], hex_vertices[:, 8]],
        [hex_vertices[:, 2], hex_vertices[:, 4], hex_vertices[:, 5], hex_vertices[:, 7]]
    ]

    volume = 0
    for tet in tets
        volume += tetrahedron_volume(tet)
    end

    return volume
end

"""
    wedge6_volume(wedge_vertices)

Calculate the volume of a wedge.

# Arguments
- `wedge_vertices`: The vertices of the wedge.
# Returns
- `volume`: The volume of the wedge.
"""
function wedge6_volume(wedge_vertices::Vector{Vector{Float64}})
    tets = [
        [wedge_vertices[1], wedge_vertices[2], wedge_vertices[3], wedge_vertices[4]],
        [wedge_vertices[2], wedge_vertices[3], wedge_vertices[4], wedge_vertices[5]],
        [wedge_vertices[3], wedge_vertices[4], wedge_vertices[5], wedge_vertices[6]]
    ]

    volume = 0
    for tet in tets
        volume += tetrahedron_volume(tet)
    end

    return volume
end

"""
    area_of_polygon(vertices)

Calculate the area of a polygon.

# Arguments
- `vertices`: The vertices of the polygon.
# Returns
- `area`: The area of the polygon.
"""
function area_of_polygon(vertices)
    n = length(vertices)
    area = 0.0

    for i in 1:n
        j = mod(i, n) + 1
        area += (vertices[i][1] + vertices[j][1]) * (vertices[i][2] - vertices[j][2])
    end

    return abs(area) / 2.0
end
