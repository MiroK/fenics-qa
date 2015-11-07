from dolfin import *
import numpy as np

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 2)
u = TrialFunction(V)
v = TestFunction(V)

facet_f = FacetFunction('size_t', mesh, 0)
AutoSubDomain(lambda x, on_boundary: x[0] > 0.5).mark(facet_f, 1)

bc = DirichletBC(V, Constant(0), facet_f, 1)

f = interpolate(Constant(1), V)
bc.apply(f.vector())

plot(f)
interactive()

# Now suppose you only want some of the facet dofs to be zeroed
dof_list = np.random.choice(np.arange(V.dim()), 20)

bc_values = bc.get_boundary_values()
for dof in dof_list:
   if dof in bc_values.keys(): del bc_values[dof]

print len(bc_values)
print len(bc.get_boundary_values())
