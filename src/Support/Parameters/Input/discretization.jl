# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct ExternalTopologyParams
    file::String = req("File"; description = "File with the external element topology")
    add_neighbor_search::Union{Nothing,Bool} = opt("Add Neighbor Search"; default = nothing)
end

@params struct SurfaceExtrusionParams
    direction::String = req("Direction"; allowed = ["X", "Y", "Z"])
    step_x::Float64 = req("Step_X"; quantity = :length)
    step_y::Float64 = req("Step_Y"; quantity = :length)
    step_z::Float64 = req("Step_Z"; quantity = :length)
    number::Int64 = req("Number"; min = 0)
end

@params struct GcodeParams
    overwrite_mesh::Bool = req("Overwrite Mesh")
    sampling::Float64 = req("Sampling"; min = 0, quantity = :length)
    width::Float64 = req("Width"; min = 0, quantity = :length)
    height::Float64 = req("Height"; min = 0, quantity = :length)
    scale::Float64 = opt("Scale"; default = 1.0)
    start_command::Union{Nothing,String} = opt("Start Command"; default = nothing)
    stop_command::Union{Nothing,String} = opt("Stop Command"; default = nothing)
    end_command::Union{Nothing,String} = opt("End Command"; default = nothing)
    blocks::Union{Nothing,Dict{String,String}} = opt("Blocks"; default = nothing)
end

@params struct BondFilterParams
    type::String = req("Type")
    normal_x::Float64 = req("Normal X")
    normal_y::Float64 = req("Normal Y")
    normal_z::Union{Nothing,Float64} = opt("Normal Z"; default = nothing)
    lower_left_corner_x::Union{Nothing,Float64} = opt("Lower Left Corner X"; default = nothing)
    lower_left_corner_y::Union{Nothing,Float64} = opt("Lower Left Corner Y"; default = nothing)
    lower_left_corner_z::Union{Nothing,Float64} = opt("Lower Left Corner Z"; default = nothing)
    bottom_unit_vector_x::Union{Nothing,Float64} = opt("Bottom Unit Vector X"; default = nothing)
    bottom_unit_vector_y::Union{Nothing,Float64} = opt("Bottom Unit Vector Y"; default = nothing)
    bottom_unit_vector_z::Union{Nothing,Float64} = opt("Bottom Unit Vector Z"; default = nothing)
    center_x::Union{Nothing,Float64} = opt("Center X"; default = nothing)
    center_y::Union{Nothing,Float64} = opt("Center Y"; default = nothing)
    center_z::Union{Nothing,Float64} = opt("Center Z"; default = nothing)
    radius::Union{Nothing,Float64} = opt("Radius"; default = nothing, min = 0)
    bottom_length::Union{Nothing,Float64} = opt("Bottom Length"; default = nothing, min = 0)
    side_length::Union{Nothing,Float64} = opt("Side Length"; default = nothing, min = 0)
    allow_contact::Bool = opt("Allow Contact"; default = false)
end

@params struct DiscretizationParams
    type::String = req("Type"; description = "Mesh format, e.g. \"Text File\" or \"Exodus\"")
    input_mesh_file::String = req("Input Mesh File")
    input_external_topology::Union{Nothing,ExternalTopologyParams} = opt("Input External Topology";
                                                                         default = nothing)
    node_sets::Dict{String,Union{Int64,String}} = opt("Node Sets";
                                                      default = Dict{String,
                                                                     Union{Int64,String}}(),
                                                      description = "Node id, list of ids, or file")
    distribution_type::Union{Nothing,String} = opt("Distribution Type"; default = nothing)
    influence_function::Union{Nothing,String} = opt("Influence Function"; default = nothing)
    surface_extrusion::Union{Nothing,SurfaceExtrusionParams} = opt("Surface Extrusion";
                                                                   default = nothing)
    bond_filters::Dict{String,BondFilterParams} = opt("Bond Filters";
                                                      default = Dict{String,BondFilterParams}())
    horizon_mesh_scaling_x::Union{Nothing,Float64} = opt("Horizon Mesh Scaling X";
                                                         default = nothing)
    horizon_mesh_scaling_y::Union{Nothing,Float64} = opt("Horizon Mesh Scaling Y";
                                                         default = nothing)
    horizon_mesh_scaling_z::Union{Nothing,Float64} = opt("Horizon Mesh Scaling Z";
                                                         default = nothing)
    gcode::Union{Nothing,GcodeParams} = opt("Gcode"; default = nothing)
end
