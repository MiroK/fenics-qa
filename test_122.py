from dolfin import *

mesh = UnitSquareMesh(10, 10)
facet_f = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], 0)').mark(facet_f, 2)

V = VectorFunctionSpace(mesh, 'CG', 1)

mu, Kpar = 1., 10.
bcs = []

# define residual 
def F(u, lmbda, v):
    I = Identity(2)
    FF = I + grad(u)
    C = FF.T*FF
    Ic = tr(C)
    JJ  = det(FF)
    psi = (mu/2.*(Ic-2 - 2*ln(JJ))+Kpar/2.*(ln(JJ))**2)*dx - lmbda*u[1]*ds(2)
    return derivative(psi, u, v)

# define variational problem
u0 = Function(V)

dE = F(u0, 0., TestFunction(V))
ddE=derivative(dE,u0,TrialFunction(V))
problem=NonlinearVariationalProblem(dE,u0,bcs,J=ddE)
solver=NonlinearVariationalSolver(problem)
