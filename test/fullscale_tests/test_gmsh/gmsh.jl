# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

import Gmsh: gmsh

# 1. Initialize gmsh
gmsh.initialize()

# 2. Define characteristic mesh size (smaller = finer mesh)
lc = 0.2

# 3. Define geometry points
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
p4 = gmsh.model.geo.addPoint(0, 1, 0, lc)

# 4. Connect points to form lines
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# 5. Create a closed line loop
line_loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])

# 6. Create a plane surface from the loop
plane_surface = gmsh.model.geo.addPlaneSurface([line_loop])

# 7. Synchronize to push geometric updates to the model
gmsh.model.geo.synchronize()

# 8. (Optional) Only mesh 2D elements (surfaces)
#gmsh.model.mesh.setDimTypes([2])

# 9. Generate the mesh
gmsh.model.mesh.generate(2)

# 10. Save to file (MSH format, viewable in gmsh, Paraview, etc.)
gmsh.write("gmsh_2d.msh")

# 7. Extrude the surface along Z axis to create a volume (a cube)
gmsh.model.geo.extrude([2, plane_surface], 0, 0, 1)

# 8. Synchronize to push geometric updates to the model
gmsh.model.geo.synchronize()

# 9. Generate the 3D mesh
gmsh.model.mesh.generate(3)

# 10. Save to file
gmsh.write("gmsh_3d.msh")

# 11. Clean up
gmsh.finalize()
