from dolfin import *
import numpy as np
import sys

field = sys.argv[1]

mesh = UnitCubeMesh(2, 2, 2)
bmesh = BoundaryMesh(mesh, 'exterior')

if field == 'scalar':
    V = FunctionSpace(mesh, 'CG', 1)
    Vb = FunctionSpace(bmesh, 'CG', 1)
    F = Expression('x[0]+x[1]', degree=1)
else:
    V = VectorFunctionSpace(mesh, 'CG', 1)
    Vb = VectorFunctionSpace(bmesh, 'CG', 1)
    F = Expression(('2*(x[0]+x[1]) + 4*(x[1]-x[2]) - x[2]',
                    'x[0]+x[1]',
                    'x[1]+x[2]'), degree=1)
fb = interpolate(F, Vb)

dofb_vb = np.array(dof_to_vertex_map(Vb), dtype=int)
size_b = fb.vector().local_size()
dofb_vb = dofb_vb[:size_b]

vb_v = np.array(bmesh.entity_map(0), dtype=int)
v_dof = np.array(vertex_to_dof_map(V), dtype=int)

f = Function(V)
array = f.vector().get_local()

# Compose the maps dof of Vb --> vertices of Vb --> vertices of V --> dof of V
in_array = fb.vector().get_local()

nsubs = V.num_sub_spaces()
if nsubs == 0:

    X = V.tabulate_dof_coordinates().reshape((-1, 3))
    Xb = Vb.tabulate_dof_coordinates().reshape((-1, 3))

    dofmap = V.dofmap()
    first, last = dofmap.ownership_range()
    mapping0 = []
    for i, xb in zip(range(size_b), Xb):
        j = np.argmin(map(lambda q: Point(*xb).distance(Point(*q)), X))
        mapping0.append(j)
        is_local = first <= dofmap.local_to_global_index(j) < last
        if not is_local:
            print 'XXXX'

    mapping = v_dof[vb_v[dofb_vb]]
else:
    # FIXME: is there a better way?
    dim = V.ufl_element().value_size()
    vb_v = np.repeat(vb_v, dim)
    maps = []
    for i in range(nsubs):
        maps.append(v_dof[i::dim][vb_v[dofb_vb[i::dim]]])
    mapping = np.array(sum(map(list, zip(*maps)), []))

# print mapping, len(array), len(mapping), size_b, mapping0
array[mapping] = in_array
f.vector().set_local(array)
f.vector().apply('insert')

print assemble(inner(f-F, f-F)*ds(domain=mesh))

foo = interpolate(F, V)
foo.vector().axpy(-1, f.vector())

print foo.vector().array()
print 'x', mapping0
