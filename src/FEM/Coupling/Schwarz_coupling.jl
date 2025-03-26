# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Schwarz_coupling
include("../Element_formulation/Lagrange_element.jl")
include("../FEM_routines.jl")
using .Lagrange_element:
                         define_lagrangian_grid_space,
                         get_recursive_lagrange_shape_functions
using .FEM_routines: get_polynomial_degree
include("../../Support/Helpers.jl")
using .Helpers:
                find_active_nodes,
                find_point_in_polygon,
                find_point_in_hexagon,
                get_ring,
                get_hexagon,
                create_centroid_and_search_radius,
                get_nearest_neighbors,
                find_point_in_element
using LinearAlgebra

export coupling_name
export init_coupling_model
export compute_coupling
function coupling_name()
    return "Schwarz"
end

function init_coupling_model(datamanager::Module, nodes, fe_params::Dict)
    @info "Coupling $(coupling_name()) is active."
    @info "Make sure the option ''Add Neighbor Search'' is true for ''Discretization''"
    @info ""

    if !haskey(fe_params["Coupling"], "Coupling Block")
        @error "In Schwarz method a block must be defined with finite elements which is also a coupling block."
        return
    end
    block_id = datamanager.get_field("Block_Id")
    fe_nodes = datamanager.get_field("FE Nodes")
    # Defines the coupling block elements as PD nodes with all features
    fe_nodes .= fe_nodes .& (block .!= fe_params["Coupling"]["Coupling Block"])
    return datamanager
end

function compute_coupling(datamanager::Module, fem_params::Dict)
    return datamanager
end

end
