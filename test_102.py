from dolfin import *
import numpy as np

# We show how integral equations can be solved with FEniCS
# We consider an integral equation of the second type and want to find 
# u: [0, 1] -> R, such that
#
# u(x) - 0.5\int_{0}^{1}(x+1)\exp{-xy}u(y)dy = f(x)   in [0, 1],
#
# where f(x) = \exp{-x} - 0.5 + 0.5\exp{-(x+1)}. The problem has an exact
# solution u(x) = exp(-x). The problem is taken from Linear Integral Equations
# by Kress (p. 158)

f = Expression('exp(-x[0])-0.5+0.5*exp(-(x[0]+1))')

def fredholm(N=32):
    'Solve the problem on mesh with N elements.'
    # Discretize the domain
    mesh = UnitIntervalMesh(N)

    # Setup Vh
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    # The discretized problem becomes M*U - 0.5*M*C*U = b
    # where M is the mass matrix, C is the matrix of the kernel, b is the load 
    # vector and U has unknown expansion coefficients of the approximate solution
    # uh

    # Mass matrix
    m_form = inner(u, v)*dx
    M = assemble(m_form).array()

    kernel = Expression('(t+1)*exp(-t*x[0])', t=0, degree=8)
    x_ks = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, 1))[:, 0]

    # Kernel matrix
    C = np.zeros(M.shape)
    # The row of C for x_k fixed can be made from the linear form
    row_form = inner(kernel, v)*dx
    for row, x_k in enumerate(x_ks):
        kernel.t = x_k
        row_values = assemble(row_form).array()
        C[row, :] = row_values

    # Load vector
    L = inner(f, v)*dx
    b = assemble(L).array()

    A = M - 0.5*M.dot(C)
    U = np.linalg.solve(A, b)

    # Create uh
    uh = Function(V)
    uh.vector()[:] = U

    # Compare with exact solution
    u_exact = Expression('exp(-x[0])')

    return errornorm(u_exact, uh, 'L2'), mesh.hmin()

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    from math import log as ln
    
    # Run the convergence test
    for N in [16, 32, 64, 128, 256, 512, 1024]:
        e, h = fredholm(N=N)
        if N != 16:
            rate = ln(e/e_)/ln(h/h_)
            print 'N=%d error=%4g rate=%.2f' % (N, e, rate)
        e_, h_ = e, h
