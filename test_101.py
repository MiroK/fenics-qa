from dolfin import *
import numpy as np

mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 10, 10)
V = FunctionSpace(mesh, 'N1curl', 1)

v = TestFunction(V)
#n = FacetNormal(mesh)

class Normal(Expression):
    def value_shape(self): return (2, )
    def eval(self, values, x):
        if near(x[0], 1):
            values[0] = 1
            values[1] = 0
        elif near(x[0], -1):
            values[0] = -1
            values[1] = 0
        elif near(x[1], 1):
            values[0] = 0
            values[1] = 1
        elif near(x[1], -1):
            values[0] = 0
            values[1] = -1
        else:
            raise ValueError("%r" % x)

n = Normal()

L = dot(v, n)*ds(metadata={'quadrature_degree': 0})
b = assemble(L)

bc = DirichletBC(V, Constant((0, 0)), 'on_boundary')
bc_dofs = sorted(bc.get_boundary_values().keys())

signs = b.array()[bc_dofs]
signs[signs > 0] = 1
signs[signs < 0] = -1
signs = dict((dof, sign) for (dof, sign) in zip(bc_dofs, signs))

for dof, sign in signs.items():
    print dof, sign
# Check

exit()

dofmap = V.dofmap()
element = V.element()
mesh.init(1)
mesh.init(1, 2)
facet_f = FacetFunction('size_t', mesh, 0)
DomainBoundary().mark(facet_f, 1)
for facet in SubsetIterator(facet_f, 1):
    cell_index = facet.entities(2)[0]
    cell = Cell(mesh, cell_index)

    dofs = dofmap.cell_dofs(cell.index())
    facet_index = cell.index(facet)
    local = dofmap.tabulate_facet_dofs(facet_index)[0]
    facet_dof = dofs[local]
    facet_mid = dofmap.tabulate_coordinates(cell)[local]
    cell_mid = cell.midpoint()
    cell_mid = np.array([cell_mid.x(), cell_mid.y()])

    vec0 = Point(facet_mid - cell_mid)

    vec1 = np.zeros(2)
    coordinate_dofs = cell.get_vertex_coordinates()
    element.evaluate_basis(local, vec1, facet_mid, coordinate_dofs,
                           cell.orientation())
    vec1 = Point(vec1)

    sign = vec1.dot(vec0)
    sign = 1 if sign > 0 else -1

    print '@', facet_dof, sign, signs[facet_dof]

