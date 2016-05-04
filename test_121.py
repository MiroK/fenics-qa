from dolfin import *

mesh = UnitSquareMesh(10, 10)
facet_f = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], 0)').mark(facet_f, 2)

ds = Measure("ds")[facet_f]

V = FunctionSpace(mesh, 'CG', 1)
c = Expression('1', domain=mesh)

# Are the measures the same?
print assemble(c*ds(2))
print assemble(c*ds(2)(mesh))

# No! So what happens
print ds(2).subdomain_id()          # 2 as expected
print ds(2)(mesh).subdomain_id()    # everywhere, why?

# Because ds(2)(mesh) is a call to ds(2).__cal__ and that results in
print ds(2).__call__.__doc__


