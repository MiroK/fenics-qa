from dolfin import *
from numpy import sqrt as nsqrt

# Define mesh and function space
mesh = IntervalMesh(500, 0, 5)
V = FunctionSpace(mesh, 'Lagrange', 1)

# Build essential boundary conditions
def u0_boundary(x, on_boundary):
    return on_boundary

# Define potential
class Potential(Expression):
    def __init__(self, a=5, b=1, c=0.5):
        Expression.__init__(self)
        self.a, self.b, self.c = a, b, c

    def eval(self, values, x):
        a, b, c = self.a, self.b, self.c
        self.M = M = 100
        # Simple quadratic potential
        #values[0] = x[0]**2
        # Double well potential 0
        values[0] = M*(abs(x-a/2.) - b)**2/b**2
        # Double well potential 1
        #values[0] = 100*((x-a/2.)**2 - b**2)**2/b**4

bc = DirichletBC(V,Constant(0.0) , u0_boundary)

# Fefine functions
u = TrialFunction(V)
v = TestFunction(V)

potential = Potential()

# Define problem
a = (inner(grad(u), grad(v)) \
     + potential*u*v)*dx
m = u*v*dx

A = PETScMatrix()
M = PETScMatrix()
_ = PETScVector()
L = Constant(0.)*v*dx

assemble_system(a, L, bc, A_tensor=A, b_tensor=_)
assemble_system(m, L, bc, A_tensor=M, b_tensor=_)

#create eigensolver
eigensolver = SLEPcEigenSolver(A,M)
eigensolver.parameters['spectrum'] = 'smallest magnitude'
eigensolver.parameters['solver'] = 'lapack'
eigensolver.parameters['tolerance'] = 1.e-15

#solve for eigenvalues
print 'Computing eigenvalues ...'
eigensolver.solve()
print 'Done.'

u = Function(V)

Au = Function(V)
Mu = Function(V)

for i in range(eigensolver.get_number_converged()):
    #extract next eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    if not near(r, 1, 1E-10):
        print i, ':', potential.M/r

        #assign eigenvector to function
        u.vector()[:] = rx #nsqrt(rx*rx)
        plot(u, interactive=True)

        Au.vector()[:] = A*u.vector()
        Mu.vector()[:] = M*u.vector()

        Au.vector().axpy(-r, Mu.vector())
        print Au.vector().max()

        rx *= -1
        print rx.sum()



