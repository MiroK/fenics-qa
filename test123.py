from dolfin import *

R = 3
h = 0.1
mesh = CircleMesh(Point(0, 0), R, h)

def top_arc(x, on_boundary):
    return on_boundary and x[1] > -DOLFIN_EPS

V = FunctionSpace(mesh, 'CG', 1)
f = Expression('0.5*(x[0]*x[0]+x[1]*x[1])', element=V.ufl_element())

facet_f = FacetFunction('size_t', mesh, 0)
AutoSubDomain(top_arc).mark(facet_f, 1)

W = VectorFunctionSpace(mesh, 'CG', 1)
gradf = project(grad(f), W)

n = FacetNormal(mesh)
ds = Measure('ds')[facet_f]
top_arc_length = assemble(dot(gradf, n)*ds(1))

print top_arc_length, pi*R*R
