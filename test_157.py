from dolfin import *
import numpy as np

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

# Suppose these will be coefficients of your form
coefs = (Expression('1', cell=mesh.ufl_cell(), degree=1),
         Expression(('x[0]', 'x[1]'), degree=1, cell=mesh.ufl_cell()),
         Constant(((1, 2), (3, 4))))

c0, c1, c2 = coefs
form = c0*inner(u, v)*dx+sqrt(c1[0]**2+c1[1]**2)*inner(u, v)*ds + det(c2)*inner(grad(u), grad(v))*dx
# All coefs in the form
print len(form.coefficients())

# Update the coefs: each coef is mapped to `ones` of the corresponding shape
update = dict((c, Constant(np.ones(c.ufl_shape))) for c in form.coefficients())
form_new = replace(form, update)
# All coefs in the form
print len(form_new.coefficients())

# Form with repeated coefficients is identfied correctly
form = c0*inner(u, v)*dx + c0*inner(u, v)*ds
assert len(form.coefficients()) == 1
