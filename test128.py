from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)

T = TrialFunction(V)
v = TestFunction(V)

F0 = (v*Dx(T,1)+Dx(v,0)*Dx(T,0))*dx - Constant(0.)*v*dx

g = Expression('x[0]+x[1]')
bc = DirichletBC(V, g, DomainBoundary())

T = interpolate(g, V)
solve(F0 == 0, T, bc,
      solver_parameters={"newton_solver":{"relative_tolerance": 1e-6,
                                          'nonzero_initial_guess': True}})

plot(T)
interactive()
