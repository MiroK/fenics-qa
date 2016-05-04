from dolfin import *
mesh = UnitSquareMesh(10, 10)
facet_f = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], 0)').mark(facet_f, 2)

ds = Measure("ds", domain=mesh, subdomain_data=facet_f)

P1 = VectorFunctionSpace(mesh, "Lagrange", 1) 
P2 = FunctionSpace(mesh, "Lagrange",1)
V = P1 * P2

mu, Kpar = 1., 10.

bc1 = DirichletBC(V.sub(0).sub(1), Constant(0.0), facet_f, 2)
bcs = [bc1]

calT = Expression(("t2"),t2=0.)

w = Function(V)
(u,p) = w.split();

I = Identity(2)
FF = I + grad(u)
C = FF.T*FF
Ic = tr(C)
JJ  = det(FF)

psi = (mu/2.*(Ic-2) - p*(JJ-1))*dx
W = psi + calT*u[0]*ds(2) + Constant(0)*inner(w, w)*dx
v = TestFunction(V)
du = TrialFunction(V)

F = derivative(W,w,v)
J = derivative(F,w,du)

problem = NonlinearVariationalProblem(F,w, bcs, J=J)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.nonlinear_solver = "snes"
solver.parameters.snes_solver.linear_solver = "umfpack"
solver.solve()
