## Dev Steps
    1. dof für koordinaten an alle cores verteilen verteilen -> done
    2. nslave nmaster verteilen -> done
    3. coordinatenfelder initialisieren auf allen Knoten -> done
    4. versenden von Koordinaten -> done
    5. versenden der Block_Id der Knoten -> done
    !-----------
    6. loc to glob -> distributionsfeld an die Knoten 
    7. glob to loc ableiten für RB Knoten
    8. overlapmap verteilen an Knoten -> send to all! -> done
    ! ---------
    9. block filter -> auf den Kernen blockID -> alle Knoten
    10. BC Filter -> an alle
    ! -------------
    11. Nachbarschaftslisten versenden ! wie in sinnvoller Weise?


## Design decisions
Each vector entry for a value exists for all nodes, also if the node does not have this property in a block. However, the synchronisation is very ugly, because all slave nodes of block with value I need the entry at the other core to. If not it will lead nowhere if MPI communication occurs

## planned feature
    search for jl files in material
    check the header for the material name
    include the file in the code via a makro
    material inclusion is very simple


