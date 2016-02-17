from dolfin import *

# create a 3*3 mesh
nx = 3
mesh = UnitSquareMesh(nx, nx)
V = FunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)

class OriginPoint(SubDomain):
    def inside(self, x, on_boundary):
        '''
        I expect to get a single vertex at the origin, but it seems not working.
        '''
        tol = 1.0/nx
        return x[0] < tol and x[1] < tol
originpoint = OriginPoint()

# create a rhs linear vector
b = assemble(v*Constant(0.0)*dx)

# create a DirichletBC
bc = DirichletBC(V, Constant(1.0), originpoint, method='pointwise')

print 'Before applying BC...'
print b.array()

bc.apply(b)
print 'After applying BC...'
print b.array()
print "changed entries:"
print bc.get_boundary_values()
