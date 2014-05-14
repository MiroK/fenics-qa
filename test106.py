from dolfin import *

set_log_level(WARNING)

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS)

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        #return x[2] < DOLFIN_EPS and on_boundary
        return near(x[2], 0) and on_boundary

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        #return x[2] > 10.0 - DOLFIN_EPS and on_boundary
        return near(x[2], 10) and on_boundary

#---------------------------------------------------------------------
def Box_chorin():

  mesh = BoxMesh(0, 0, 0, 1, 1, 10, 12, 12, 60)
  h = mesh.hmin()

  #! problem specific
  f = Constant((0., 0., 0.))     # force
  Se = Constant(0.)              # energy source
  nu = Constant(1./8.)           # kinematic viscosity
  lamda = Constant(0.6)          # Thermal Conductivity
  Cp = Constant(4200.)           # Specific Heat Capacity
  dt = 0.2*h/1.                  # time step CFL with 1 = max. velocity
  k = Constant(dt)               # time step
  Time = 1.0                     # total simulation time
  u0 = Constant((0., 0., 0.))    # initial velocity
  p0 = Constant(0.)              # initial pressure
  T0 = Constant(303.)            # initial temperature = 303 K

  #! solver specific
  V = VectorFunctionSpace(mesh, 'CG', 2)
  Q = FunctionSpace(mesh, 'CG', 1)


  u = TrialFunction(V)
  v = TestFunction(V)
  p = TrialFunction(Q)
  q = TestFunction(Q)
  T = TrialFunction(Q)
  Tq = TestFunction(Q)

  u0 = interpolate(u0, V)
  p0 = interpolate(p0, Q)
  T0 = interpolate(T0, Q)
  us = Function(V)
  u1 = Function(V)
  p1 = Function(Q)
  T1 = Function(Q)

  # tentative velocity
  F0 = (1./k)*inner(u - u0, v)*dx + inner(dot(grad(u0), u0), v)*dx\
       + nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
  a0, L0 = system(F0)

  # projection
  F1 = inner(grad(p), grad(q))*dx + (1./k)*q*div(us)*dx
  a1, L1 = system(F1)

  # finalize
  F2 = (1./k)*inner(u - us, v)*dx + inner(grad(p1), v)*dx
  a2, L2 = system(F2)

  # Energy
  F3 = Cp*(1./k)*inner(T - T0, Tq)*dx + Cp*inner(dot(grad(T0), u1), Tq)*dx\
       + lamda*inner(grad(T), grad(Tq))*dx - inner(Se, Tq)*dx
  a3, L3 = system(F3)

  # boundary conditions
  # Walls
  b_v = DirichletBC(V, Constant((0.0, 0.0, 0.0)), NoslipBoundary())
  b_T = DirichletBC(Q, Constant(323.0), NoslipBoundary())       # Walls T = 323 K

  # Inlet
  b_v1 = DirichletBC(V, Constant((0.0, 0.0, 0.3)), InflowBoundary())
  b_T1 = DirichletBC(Q, Constant(303.0), InflowBoundary())      # Inlet T = 303 K

  # Outlet
  b_p1 = DirichletBC(Q, Constant(0.), OutflowBoundary())

  bcs_v = [b_v, b_v1]
  bcs_p = [b_p1]
  bcs_T = [b_T, b_T1]

  A0 = assemble(a0)
  A1 = assemble(a1)
  A2 = assemble(a2)
  A3 = assemble(a3)

  solver02 = KrylovSolver('gmres', 'ilu')
  solver1 = KrylovSolver('cg', 'petsc_amg')
  solver3 = KrylovSolver('gmres', 'ilu')

  ufile = File("Box_chorin_u.pvd")
  pfile = File("Box_chorin_p.pvd")
  Tfile = File("Box_chorin_T.pvd")

  iter = 0
  t = 0
  while t < 2*dt:
    t += dt
    iter += 1

    b = assemble(L0)
    [bc.apply(A0, b) for bc in bcs_v]
    solver02.solve(A0, us.vector(), b)

    b = assemble(L1)
    [bc.apply(A1, b) for bc in bcs_p]
    solver1.solve(A1, p1.vector(), b)

    b = assemble(L2)
    [bc.apply(A2, b) for bc in bcs_v]
    solver02.solve(A2, u1.vector(), b)

    b = assemble(L3)
    [bc.apply(A3, b) for bc in bcs_T]
    solver3.solve(A3, T1.vector(), b)

    ufile << u1
    pfile << p1
    Tfile << T1

    u0.assign(u1)
    T0.assign(T1)
    plot(T0, interactive=True)
    print t
    print iter

  return u0, p1, T0
#----------------------------------------------------------------------
if __name__ == '__main__':
  u, p, T = Box_chorin()
