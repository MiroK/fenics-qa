from dolfin import *

mesh = UnitCubeMesh(40, 40, 40)
tdim = mesh.topology().dim()

mesh.init(tdim-1, 0)
v_2_f = dict((tuple(facet.entities(0)), facet.index())
             for facet in facets(mesh))

print len(v_2_f)
