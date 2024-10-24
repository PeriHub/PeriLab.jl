# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Arlequin_coupling
function coupling_name()
    return "Arlequin"
end

function init_coupling_model(datamanager::Module, complete_params::Dict)

    datamanager.create_constant_node_field("Number of Element Neighbors", Int64, 1, 1)
    # Find number of element neighbors
    datamanager.create_constant_element_field("Coupling Elementlist", Int64, 1)
    # Assign Element ids

    return datamanager
end

function compute_coupling(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    fem_params::Dict,
)

    coupling_elements = datamanager.get_field("Coupling Elementlist") # elements connected to PD nodes
    # p1->(E1,E4)->E1-> 1,3,7,2..
    topology = datamanager.get_field("FE Topology")
    coordinates = datamanager.get_field("Coordinates")
    displacements = datamanager.get_field("Displacements", "NP1")
    #lumped_mass => hier rechnen oder nicht
    # Jinv = datamanager.get_field("Inverse Jacobian")
    # dichte oder Massenmatrix, tbthink->als funktion
    N_Matrix = datamanager.get_field("N Matrix") #-> mapping im element; eventuell NT(NTN)^-1 speichern
    # das optional
    B_Matrix = datamanager.get_field("B Matrix")
    # position des PD Knoten im element wird hier gemacht
    #belieige postionen im Element

    #distribut_forces(dasd)
    @info nodes
    @info coupling_elements
    @info topology
    @info coordinates
    @info displacements
    # @info Jinv
    @info N_Matrix
    @info B_Matrix
    return datamanager

end




end
