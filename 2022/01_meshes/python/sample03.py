import gmsh
import sys

gmsh.initialize()

gmsh.model.add("t3")

lc = 1e-2
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(.1, 0, 0, lc, 2)
gmsh.model.geo.addPoint(.1, .3, 0, lc, 3)
gmsh.model.geo.addPoint(0, .3, 0, lc, 4)

gmsh.model.geo.addPoint(.02, .02, 0, lc, 5)
gmsh.model.geo.addPoint(.08, 0.02, 0, lc, 6)
gmsh.model.geo.addPoint(.08, .28, 0, lc, 7)
gmsh.model.geo.addPoint(0.02, .28, 0, lc, 8)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(3, 2, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

gmsh.model.geo.addLine(5, 6, 5)
gmsh.model.geo.addLine(7, 6, 6)
gmsh.model.geo.addLine(7, 8, 7)
gmsh.model.geo.addLine(8, 5, 8)

gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)
gmsh.model.geo.addCurveLoop([8, 5, -6, 7], 2)

gmsh.model.geo.addPlaneSurface([1, -2], 1)

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("t3.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

