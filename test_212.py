from dolfin import *
import numpy

h = []  # mesh sizes
E = []  # errors

f=Expression('sin(pi*x[1])*(pi*pi*x[0])', degree = 10)

for nx in [4,8,16,32,64]:
    h.append(1.0/nx)
    mesh = UnitSquareMesh(nx, nx)
    Vc = FiniteElement("CG", mesh.ufl_cell(), 2)
    V0 = FiniteElement("DG", mesh.ufl_cell(), 0)
    Vh = FunctionSpace(mesh, MixedElement([Vc, V0]))

    u_exact=Expression('x[0]*sin(pi*x[1])', degree = 6, domain = mesh)

    u, u0 = TrialFunction(Vh)
    v, v0 = TestFunction(Vh)
    n = FacetNormal(mesh)
    hcell = CellSize(mesh)
    a = inner(grad(u), grad(v))*dx \
    + (- inner(avg(grad(u)), n("+"))*jump(v0) - inner(avg(grad(v)), n("+"))*jump(u0) \
    +  (10.0/avg(hcell))*jump(u0)*jump(v0) )*dS \
    + (- inner(grad(u), n)*(v+v0) - inner(grad(v), n)*(u+u0))*ds + (10.0/hcell)*(u+u0)*(v+v0)*ds
    b = inner(f, v)*dx + inner(f, v0)*dx \
    + 10./hcell*u_exact*(v+v0)*ds - u_exact*dot(grad(v), n)*ds

    u_sol = Function(Vh)
    solve(a == b, u_sol)
    (u, u0) = u_sol.split(deepcopy=True)

    # Compute errors        
    uH1_error= errornorm(u_exact, u, norm_type = 'H10', degree_rise=3)
    uH1man_error= sqrt(assemble(dot(grad(u_exact) - grad(u), grad(u_exact) - grad(u))*dx))
    errors = {'u_H1_error' : uH1_error, 'u_H1man_error' : uH1man_error}
    E.append(errors)

# Print convergence rates
from math import log as ln
error_types = E[0].keys()
for error_type in sorted(error_types):
    j = len(E)
    print '\nError norm of', error_type
    for i in range(0, j-1):
        r = ln(E[i+1][error_type]/E[i][error_type])/ln(0.5)  # E is a list of errors
        print 'mesh size =%.7f Ei=%.7f Ei+1=%.7f r=%.2f' % (h[i], E[i][error_type], E[i+1][error_type], r)
