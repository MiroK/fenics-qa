from dolfin import *
import numpy as np

manual = True  # How to set bcs

assert MPI.size(mpi_comm_world()) == 1, 'The code below is intended to run in serial'

mesh = UnitSquareMesh(40, 40)

# Build the set of vertices where bcs should be prescribed
node_set = VertexFunction('size_t', mesh, 0)
bc_boundary = CompiledSubDomain('near(x[0], 0)')
bc_boundary.mark(node_set, 1)
# Only keep vertices marked as 1
node_set = [v.index() for v in SubsetIterator(node_set, 1)]

# Get dofs corresponsing to vertices
V = FunctionSpace(mesh, 'CG', 1)
dof_set = np.array(vertex_to_dof_map(V)[node_set], dtype='intc')

# Assemble the system
u, v = TrialFunction(V), TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = inner(Constant(1), v)*dx
A, b = assemble_system(a, L)

bc_f = Constant(2)  # Boundary value to be prescribed
# Manual application of bcs
if manual:
    # Modif A: zero bc row & set diagonal to 1
    A.ident_local(dof_set)
    A.apply('insert')

    # Modif b: entry in the bc row is taken from bc_f
    bc_values = interpolate(bc_f, V).vector().array()
    b_values = b.array()
    b_values[dof_set] = bc_values[dof_set]
    b.set_local(b_values)
    b.apply('insert')
# Auto at same domain which determined node_set
else:
    bc = DirichletBC(V, bc_f, 'near(x[0], 0)')
    bc.apply(A)
    bc.apply(b)

# Check that auto matches manual
uh = Function(V)
solve(A, uh.vector(), b)
plot(uh)
interactive()

print uh.vector().norm('l2')
