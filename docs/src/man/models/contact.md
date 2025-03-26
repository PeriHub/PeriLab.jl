## Contact
## Contact Search

**Initialization**

Step 1
Identification of all surface nodes of contact blocks.

Step 2
All surface nodes are stored in a list. This list exists at all cores with the current position of the surface nodes.

Step 3
Connect each surface node to a geometrical surface.

---

**Computation**

Step 1
Perform a nearest neighbor search with a user defined search radius. The result is a list of potential contact pairs.

Step 2
The surfaces and the polyeder is updated, due to deformation and the surface normals can be computed.

Step 3
Check if the master point lies inside the polyeder. If so connect this point with the surface normal and its nearest neighbor.

----

## Contact
A bond-based contact formulation is computed. The stiffness is defined in YAML input file. The distance between master and slave surface is used. The resulting force is applied to the slave node and it's surface neighbors. The sum is equal to the force applied to the master node.
