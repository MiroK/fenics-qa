from dolfin import *

mesh = UnitSquareMesh(20, 20)
h = mesh.hmin()

sub_cells = CellFunction('size_t', mesh, 0)
AutoSubDomain(lambda x: x[0] < 0.5 + h).mark(sub_cells, 1)

sub_facets = FacetFunction('size_t', mesh, 0)
mesh.init(2, 1)
for cell in SubsetIterator(sub_cells, 1):
    for facet in facets(cell):
        sub_facets[facet] = 1

V = FunctionSpace(mesh, 'CG', 1)
f = interpolate(Expression('x[0]*x[0]+x[1]*x[1]+1'), V)
bc = DirichletBC(V, Constant(1), sub_facets, 1)
bc.apply(f.vector())

plot(f)
interactive()
