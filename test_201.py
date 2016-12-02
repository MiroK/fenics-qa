from dolfin import *
import numpy as np

npoints = 10

comm = mpi_comm_world().tompi4py()

values = np.zeros(npoints)
points = np.zeros((npoints, 2))
if comm.rank == 0:
    points = np.random.rand(npoints, 2)
    values = np.random.rand(npoints)
values = comm.bcast(values)
points = comm.bcast(values)
points = map(Point, points)

# Sources
f = Constant(0)
g = Expression(tuple(map(str, values)), degree=1)

def delta_h(point):
    '''A dirac like function localized at point'''
    return Expression('(sqrt(pow((x[0]-Px), 2)+pow((x[1]-Py), 2)) < h) ? 1 : 0',
                      degree=2, h=mesh.hmin(), Px=point.x(), Py=point.y())

ph0 = None
for ncells in [8, 16, 32, 64, 128, 256, 512]:
    mesh = UnitSquareMesh(ncells, ncells)

    # Variational form
    V = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
    Q = VectorElement('Real', mesh.ufl_cell(), 0, npoints)
    W = MixedElement([V, Q])

    W = FunctionSpace(mesh, W)

    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)

    h = CellVolume(mesh)
    a = inner(grad(u), grad(v))*dx +\
        sum((1./h)*inner(delta_h(xi)*pi, v)*dx for xi, pi in zip(points, p)) +\
        sum((1./h)*inner(delta_h(xi)*qi, u)*dx for xi, qi in zip(points, q))
    L = inner(f, v)*dx + inner(g, q)*dx

    bc = DirichletBC(W.sub(0), Constant(0), 'on_boundary')

    wh = Function(W)
    solve(a == L, wh, bc)

    uh, ph = wh.split(deepcopy=True)
    
    ph = ph.vector()
    if ph0 is not None:
        e = (ph - ph0).norm('l2')/ph0.norm('l2')
        print e
    ph0 = ph

plot(uh, title='ncells %d' % ncells)
interactive()
