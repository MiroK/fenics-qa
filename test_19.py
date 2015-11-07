from dolfin import *

parameters['linear_algebra_backend'] = 'PETSc'

mesh = UnitSquareMesh(4, 4)
V = FunctionSpace(mesh, 'CG', 1)

u1 = interpolate(Expression('1+x[0]+x[1]'), V)
u2 = interpolate(Expression('1+x[0]'), V)
uerr = project(u1 - u2,V)    
uerr_rel = project(uerr/u1, V)

print uerr_rel.vector().array()

Uerr = Function(V)
Uerr.assign(FunctionAXPY([(1., u1), (-1., u2)]))

Uerr_rel = Function(V)
Uerr_rel.vector().axpy(1, Uerr.vector())
as_backend_type(Uerr_rel.vector()).vec().pointwiseDivide(
        as_backend_type(Uerr_rel.vector()).vec(),
        as_backend_type(u1.vector()).vec())

print Uerr_rel.vector().array()


