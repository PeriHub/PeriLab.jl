<!--
SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

## Dev Steps
  Introduction of bond associated correspondence 

**Developement plan**
- Bond associated neighborhood is the overlap between nlist[iID] and nlist[nlist[iID][jID]]
- Filter equal nodes and create a new neighborhoodlist for bond -> bond_nlist
- calculate K, Kinv and defGrad -> already there if the neighborhood loop is in a function
- weighted volume (sum(volume(bond_nlist))/sum(volume[nlist[iID]]))

If this works for one core the following will be introduced

- all neighbors search for neighbors at each core
- numbers are correct and it allows a change in size -> local ID is correct


## Design decisions
Each vector entry for a value exists for all nodes, also if the node does not have this property in a block. However, the synchronisation is very ugly, because all responder nodes of block with value I need the entry at the other core to. If not it will lead nowhere if MPI communication occurs

    IO
    nodesets are not defined yet in Exodus.jl

    snake_case for variables and functions
    PascalCase for modules and type names
    FULL_UPPERCASE for constants

## Issues
    for n=4 -> errors    
    MPI_Neighbor_alltoall -> might be easier
## planned feature
    test if blocks are defined in yaml, but missing in mesh
    https://github.com/StephenVavasis/Unroll.jl
    static arrays package -> speed up -> only for arrays smaller 100; inverte of Jacobian, etc.
    integration of FEM Julia package -> coupling might be better, becaus of more functionality
    time step minimum for all cores -> parallel computing! -> done
    matrix -> reshape from vector for better use
    search for jl files in material
    check the header for the material name
    include the file in the code via a makro
    material inclusion is very simple
    bonds as elements in exodus -> filter small to large writing ?! -> elements can be x,y to be represented both
    multiple materials in one block -> evaluation order



