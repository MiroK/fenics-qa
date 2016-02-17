from dolfin import *

mesh = UnitSquareMesh(20, 20)

V = FunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)

f = Expression('x[1]')
facet_f = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], 0)').mark(facet_f, 1)
ds = Measure('ds', domain=mesh, subdomain_id=1, subdomain_data=facet_f)

whole = assemble(f*f*ds)
pieces = assemble(f*v*ds)

print whole - pieces.sum()
