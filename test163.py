from dolfin import *
from random import uniform
from numpy import zeros

comm = mpi_comm_world()

class InitialCondition(Expression):
    def eval(self, values, x):
        values[0] = .5+1*uniform(-.1,.1)

mesh = RectangleMesh(-25,-25,25,25, 25, 25, 'crossed')

V = FunctionSpace(mesh, "Lagrange", 1)
V_vec = VectorFunctionSpace(mesh, "Lagrange", 1)

u = TrialFunction(V_vec)
test_u = TestFunction(V_vec)
u_old = Function(V_vec)

def eps(u):
            return as_vector([u[i].dx(i) for i in range(2)] +
                     [u[i].dx(j) + u[j].dx(i) for (i,j) in [(0,1)]])

top = CompiledSubDomain("(std::abs(x[1]-d) < DOLFIN_EPS) && on_boundary", d = 25)


boundaries = MeshFunction('size_t', mesh, 1) # Gets object of boundaries
boundaries.set_all(0)
top.mark(boundaries, 1)
ds = Measure('ds')[boundaries]

au = inner(eps(u), eps(test_u))*dx
Lu = Constant(.1)*test_u[0]*ds(1)
u = Function(V_vec)

A, b = assemble_system(au, Lu)

file_u = XDMFFile(comm, 'output/disp.xdmf')
file_u << (u, 0.)

solver = PETScKrylovSolver('gmres')

x = mesh.coordinates()
xp = zeros(x.shape)
xp[:,0] = x[:,1]
xp[:,1] = -x[:,0]
d = x-xp

null_vec = Vector(u.vector())
null_vec2 = Vector(u.vector())
V_vec.sub(0).dofmap().set(null_vec, 1)
V_vec.sub(1).dofmap().set(null_vec2, 1)

r = Function(V_vec)
rx = Function(V)
ry = Function(V)
rx.vector().set_local(d[:,0])
ry.vector().set_local(d[:,1])
assign(r.sub(0), rx)
assign(r.sub(1), ry)
null_vec3 = r.vector()

null_vec *= 1.0/null_vec.norm("l2")
null_vec2 *= 1.0/null_vec2.norm("l2")
null_vec3 *= 1.0/null_vec3.norm("l2")

ns = [null_vec, null_vec2, null_vec3]

n0 = Constant((1, 0))
n1 = Constant((0, 1))
n2 = Expression(('x[0]-x[1]', 'x[1]-x[0]'))
ns = [n0, n1, n2]

# Basis of rigid motions in V_vec
null_space_basis = [interpolate(n, V_vec).vector() for n in ns]

# Make into unit vectors
[normalize(n, 'l2') for n in null_space_basis]
ns = null_space_basis
for i, ei in enumerate(ns):
    for j, ej in enumerate(ns):
        print i, j, ei.inner(ej)



# Create null space basis object and attach to Krylov solver
null_space = VectorSpaceBasis(null_space_basis)
solver.set_nullspace(null_space)

 # Orthogonalize b with respect to the null space (this gurantees that
 # a solution exists)
null_space.orthogonalize(b);

solver.solve(A, u.vector(), b)
file_u << (u, 1.)

plot(u, interactive = True)
