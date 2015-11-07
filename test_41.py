from dolfin import *

nx = 64
ny = 64
mesh = UnitSquareMesh(nx, ny)

V = FunctionSpace(mesh, "CG", 1) 
R = FunctionSpace(mesh, "R", 0)
W = V * R
(u, c) = TrialFunction(W)
(v, d) = TestFunction(W)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
g = Expression("-sin(5*x[0])")
a = ( inner(grad(u),grad(v)) + c*v + d*u )*dx
L = f*v*dx + g*v*ds

w = Function(W)
n_prob_laplace_operator = LinearVariationalProblem( a, L, w )
solver_np_laplace_operator = LinearVariationalSolver(n_prob_laplace_operator)
solver_np_laplace_operator.solve()
(u, c) = w.split()


# BEGIN: PROBLEMATIC BLOCK OF CODE
# Dump solution to the screen
u_nodal_values = u.vector()
u_array = u_nodal_values.array()
coor = mesh.coordinates()
if coor.shape[0] == u_array.shape[0]:  # degree 1 element
    for i in range(len(u_array)):
        print 'u(%8g,%8g) = %g' % (coor[i][0], coor[i][1], u_array[i])
else:
    print """\
 Cannot print vertex coordinates and corresponding values
 because the element is not a first-order Lagrange element.
 """

# END: PROBLEMATIC BLOCK OF CODE

plot(u, interactive=True)
