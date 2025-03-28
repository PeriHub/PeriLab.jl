## Contact
## Contact Search

| Parameter | Unit | Description |
|---|---|---|
|Search Radius | $[m]$| Radius to get a list of potential contact pairs.|
|Master| $[-]$| Block ID of the master nodes. |
|Slave | $[-]$| Block ID of the slave nodes.|
**Initialization**

- Step 1 -
Identification of all surface nodes of contact blocks.
- Step 2 -
All surface nodes are stored in a list. This list exists at all cores with the current position of the surface nodes.
- Step 3 -
Connect each surface node to a geometrical surface.

**List mapping**
To fill the postion vector with the needed values mappings are needed.

local point ids -> global point id -> reduced contact block point id
The mapping has to be done in both directions, because the contact forces have to be applied to the local core point id.



---

**Computation**

- Step 1-
Perform a nearest neighbor search with a user defined search radius. The result is a list of potential contact pairs.
- Step 2 -
The surfaces and the polyeder is updated, due to deformation and the surface normals can be computed.
- Step 3 -
Check if the master point lies inside the polyeder. If so connect this point with the surface normal and its nearest neighbor.

**Synchronization**
All cores (except core 1) send it's displacement values to core 1. Core 1 sends them back.
!!! note "Efficiency"
    Smaller contact areas are more efficient in numerical analysis, because less synchronisation have to occur.

----

## Contact
A bond-based contact formulation is computed. The stiffness is defined in YAML input file. The distance between master and slave surface is used. The resulting force is applied to the slave node and it's surface neighbors. The sum is equal to the force applied to the master node.
