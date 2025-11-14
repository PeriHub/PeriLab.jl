# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
start julia
]
activate .
run the code in the main structure in the REPL
if not in the REPL, make sure you that you adapt the pathes

"""
# for REPL
include("./src/Core/Data_manager.jl")

using .Data_manager

dm = Data_manager
dm.initialize_data()

# 3x3 Grid Setup (9 nodes total)
nodes = 9
dof = 2  # 2D problem

# Initialize data manager
dm.initialize_data()
dm.set_num_controller(nodes)
dm.set_dof(dof)

# Create node fields
nn = dm.create_constant_node_field("Number of Neighbors", Int64, 1)
h = dm.create_constant_node_field("Horizon", Float64, 1)
volume = dm.create_constant_node_field("Volume", Float64, 1, 0.1)
dm.create_node_field("Forces", Float64, dof)
dm.create_node_field("Force Densities", Float64, dof)
volume = dm.create_constant_node_field("Force", Float64, 1, 0.1)

# Define 3x3 grid positions (spacing = 1.0)
# Node layout:
# 7 8 9
# 4 5 6
# 1 2 3

positions = [
    [0.0, 0.0],  # Node 1
    [1.0, 0.0],  # Node 2
    [2.0, 0.0],  # Node 3
    [0.0, 1.0],  # Node 4
    [1.0, 1.0],  # Node 5 (center)
    [2.0, 1.0],  # Node 6
    [0.0, 2.0],  # Node 7
    [1.0, 2.0],  # Node 8
    [2.0, 2.0]   # Node 9
]
dg_N, dg_NP1 = dm.create_node_field("Deformed Coordinates", Float64, dof)
dg_N = copy(positions)
# Set horizon for each node (e.g., 1.5 to reach nearest neighbors)
h[1:nodes] .= 1.5

# Define number of neighbors for each node
# Corner nodes: 2 neighbors, Edge nodes: 3 neighbors, Center node: 4 neighbors
nn[1] = 2  # Corner
nn[2] = 3  # Edge
nn[3] = 2  # Corner
nn[4] = 3  # Edge
nn[5] = 4  # Center
nn[6] = 3  # Edge
nn[7] = 2  # Corner
nn[8] = 3  # Edge
nn[9] = 2  # Corner

# Create bond fields
bf = dm.create_constant_bond_field("Bond Forces", Float64, dof)
bdN, bdNP1 = dm.create_bond_field("Bond Damage", Float64, 1, 1)
dbN, dbNP1 = dm.create_bond_field("Deformed Bond Geometry", Float64, dof, 1)
dbdN, dbdNP1 = dm.create_bond_field("Deformed Bond Length", Float64, 1)
bg = dm.create_constant_bond_field("Bond Geometry", Float64, dof)
bd = dm.create_constant_bond_field("Bond Length", Float64, 1, 1)
nlist = dm.create_constant_bond_field("Neighborhoodlist", Float64, 1, 1)

nlist = [
    [2, 4],        # Node 1: neighbors 2, 4
    [1, 3, 5],     # Node 2: neighbors 1, 3, 5
    [2, 6],        # Node 3: neighbors 2, 6
    [1, 5, 7],     # Node 4: neighbors 1, 5, 7
    [2, 4, 6, 8],  # Node 5: neighbors 2, 4, 6, 8 (center)
    [3, 5, 9],     # Node 6: neighbors 3, 5, 9
    [4, 8],        # Node 7: neighbors 4, 8
    [5, 7, 9],     # Node 8: neighbors 5, 7, 9
    [6, 8]         # Node 9: neighbors 6, 8
]
