from dolfin import *

def eps(u):
    return as_vector([u[i].dx(i) for i in range(2)] +
                     [u[i].dx(j) + u[j].dx(i) for (i,j) in [(0,1)]])

fcn = Expression('x[0]*.01')

mesh = RectangleMesh(0,0,1,1, 10, 10) #, 'crossed')
V1 = FunctionSpace(mesh, "Lagrange", 1)
V_vec = VectorFunctionSpace(mesh, "Lagrange", 1)

c = Function(V1)

u = TrialFunction(V_vec)
test_u = TestFunction(V_vec)

left = CompiledSubDomain("(std::abs(x[0]) < DOLFIN_EPS) && on_boundary")

bc = DirichletBC(V_vec, [0, 0] , left)
eps0 = as_vector((1,1,0))

c.interpolate(fcn)

au = inner(eps(u), eps(test_u))*dx
Lu = c*inner(eps0, eps(test_u))*dx

B = c*inner(u, eps(test_u))*dx

u = Function(V_vec)

m_el = inner(eps0, eps(test_u))*dx
M_el = assemble(m_el)

A = assemble(au)
bc.apply(A)

x = B*project(as_vector((c,c)),V_vec).vector()
print x.norm()
b = assemble(Lu) # Uncomment to see the correct answer
print b.norm()
bc.apply(b)


solver_u = PETScKrylovSolver('cg')
prm = solver_u.parameters
prm["absolute_tolerance"] = 1e-4
prm["relative_tolerance"] = 1e-5
prm["nonzero_initial_guess"] = True

solver_u.solve(A, u.vector(), b)

plot(u, interactive=True)
