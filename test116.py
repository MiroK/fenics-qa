from dolfin import *

def save_gnuplot(file_name, u):
    mesh = u.function_space().mesh()
    t_dim = mesh.topology().dim()
    g_dim = mesh.geometry().dim()
    vertex_coordinates = numpy.zeros((t_dim+1)*gdim)
    
    with open(file_name, 'w') as gnuplot_file:
        for i, cell in enumerate(cells(mesh)):
            cell.get_vertex_coordinates(vertex_coordinates)
            for 

# Export to gnuplot
def gnuplot_export(file, mesh, function):
  file = open(file, 'w+')
  i = 0
  for myCell in cells(mesh):
    i += 1
    myVertices = vertices(myCell)
    it = iter(myVertices)
    myVertex0 = it.next()
    myVertex1 = it.next()
    myVertex2 = it.next()
    print >>file, '%e  %e  %e' % (myVertex0.x(0),myVertex0.x(1),function(myVertex0.x(0),myVertex0.x(1)))
    print >>file, '%e  %e  %e' % (myVertex1.x(0),myVertex1.x(1),function(myVertex1.x(0),myVertex1.x(1)))
    print >>file, '%e  %e  %e' % (myVertex2.x(0),myVertex2.x(1),function(myVertex2.x(0),myVertex2.x(1)))
    print >>file, '%e  %e  %e' % (myVertex0.x(0),myVertex0.x(1),function(myVertex0.x(0),myVertex0.x(1)))
    if (i == 1):
          print >>file, '%e  %e  %e' % (myVertex0.x(0),myVertex0.x(1),function(myVertex0.x(0),myVertex0.x(1)))
    print >>file, ''

NUM_CELL = 24
mesh = UnitSquareMesh(NUM_CELL,NUM_CELL)
h = CellSize(mesh)
n = FacetNormal(mesh)

# Create FunctionSpaces
V = FunctionSpace(mesh, "CG", 3)
W = FunctionSpace(mesh, "CG", 1)

# Boundary conditions
def right(x, on_boundary): return x[0] > (1.0 - DOLFIN_EPS)
def left(x, on_boundary): return x[0] < DOLFIN_EPS
def bottom_center(x, on_boundary):
    return x[1] < DOLFIN_EPS and (x[0] > 1./3. - DOLFIN_EPS and x[0] < 2./3.0 + DOLFIN_EPS)
def bottom_lr(x, on_boundary):
    return x[1] < DOLFIN_EPS and (x[0] < 1./3. + DOLFIN_EPS or x[0] > 2./3.0 - DOLFIN_EPS)
def top(x, on_boundary):
    return x[1] > 1.0 - DOLFIN_EPS

g0 = Constant(0.0)
g1 = Constant(1.0)

bc0 = DirichletBC(V, g1, bottom_center)
bc1 = DirichletBC(V, g0, bottom_lr)
bc3 = DirichletBC(V, g0, top)
bc4 = DirichletBC(V, g0, right)
bcs = [bc0, bc1, bc4, bc3]

# Parameters
epsilon = Constant(0.000000001)
c = Constant(0.)
b = Expression(('-x[1]', 'x[0]'))
f = Constant(0.)

uh = TrialFunction(V)
vh = TestFunction(V)
wh = TestFunction(W)

bb = assemble(dot(b,b)*wh*dx)
bF = Function(W,bb)
tau = 1.

# SUPG (SDFEM) method
a1 = (epsilon*dot(grad(uh),grad(vh)) + vh*dot(b,grad(uh)) + c*uh*vh)*dx
a2 = (h/(2.*sqrt(bF))*tau*inner(dot(b,grad(uh)),dot(b,grad(vh))))*dx
a = a1 + a2
L = f*vh*dx + h/(2.*sqrt(bF))*tau*h*f*dot(b, grad(vh))*dx

# Compute solution
uh = Function(V)
solve(a == L, uh, bcs)

gnuplot_export('./testfile', mesh, uh)
