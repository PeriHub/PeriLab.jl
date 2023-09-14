## Dev Steps
    1. dof für koordinaten an alle cores verteilen verteilen -> done
    2. nslave nmaster verteilen -> done
    3. coordinatenfelder initialisieren auf allen Knoten -> done
    4. versenden von Koordinaten -> done
    5. versenden der Block_Id der Knoten -> done
    5.a neighborhoodlist -> num neighbor für jeden knoten -> done
                         -> feld init an jedem knoten -> das vielleicht skalierbare referenz? -> done
    5.b neighboorhoodlist verteilen -> done
    5.c bond init lists -> done
    !-----------
    6. loc to glob -> distributionsfeld an die Knoten -> done 
    7. glob to loc ableiten für RB Knoten -> done
    8. overlapmap verteilen an Knoten -> send to all! -> done
    ! ---------
    9. block filter -> auf den Kernen blockID -> alle Knoten -> done
    10. BC Filter -> an alle -> done
    ! -------------
    11. Nachbarschaftslisten versenden ! wie in sinnvoller Weise?
    12. bondvectors set to zero -> done
    13. bc interpreter -> done
    13.a nodesets in mesh -> done
    14. overlap synchronisation
    15. bc in solver -> done
    16. verlet solver -> pacackage?
    17. write output -> done
    18. integrate first model
    18.a step width determination -> done
    19. first test
    20. 2D arrays in fields
    21. params reader for material + physics -> done
    22. compute class

## Design decisions
Each vector entry for a value exists for all nodes, also if the node does not have this property in a block. However, the synchronisation is very ugly, because all slave nodes of block with value I need the entry at the other core to. If not it will lead nowhere if MPI communication occurs

    IO
    nodesets are not defined yet in Exodus.jl

## Issues
    for n=4 -> errors    
## planned feature
    test if blocks are defined in yaml, but missing in mesh
    https://github.com/StephenVavasis/Unroll.jl
    static arrays package -> speed up
    integration of FEM Julia package -> coupling might be better, becaus of more functionality
    time step minimum for all cores -> parallel computing! -> done
    matrix -> reshape from vector for better use
    search for jl files in material
    check the header for the material name
    include the file in the code via a makro
    material inclusion is very simple
    bonds as elements in exodus -> filter small to large writing ?! -> elements can be x,y to be represented both
    multiple materials in one block -> evaluation order



