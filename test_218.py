from dolfin import *
import numpy as np

mesh = UnitCubeMesh(1, 1, 1)
bmesh = BoundaryMesh(mesh, 'exterior')

V = VectorFunctionSpace(mesh, 'CG', 1)
Vb = VectorFunctionSpace(bmesh, 'CG', 1)

#F = Expression(('2*(x[0]+x[1]) + 4*(x[1]-x[2]) - x[2]',
# 3                'x[0]+x[1]',
#                'x[1]+x[2]'), degree=1)
#fb = interpolate(F, Vb)

Xb = Vb.tabulate_dof_coordinates().reshape((-1, 3))

dofb_vb = np.array(dof_to_vertex_map(Vb), dtype=int)

vb_v = np.repeat(np.array(bmesh.entity_map(0), dtype=int), 3)
print vb_v
print np.array(bmesh.entity_map(0))

v_dof = np.array(vertex_to_dof_map(V), dtype=int)

print v_dof[vb_v[dofb_vb]]
X = V.tabulate_dof_coordinates().reshape((-1, 3))
print X
#f = Function(V)
#array = f.vector().array()
# Compose the maps dof of Vb --> vertices of Vb --> vertices of V --> dof of V
#in_array = fb.vector().array()
#in_array = in_array[np.r_[Vb.sub(0).dofmap().dofs(),
#                          Vb.sub(1).dofmap().dofs(), 
#                          Vb.sub(2).dofmap().dofs()]]


for i, j in zip(X[v_dof[vb_v[dofb_vb]]], Xb):
    print 'x', i, j
#f.vector()[:] = array

#print assemble(inner(f-F, f-F)*ds(domain=mesh))

# plot(f)
# interactive()
# FIgure this thing for vector and parallel
