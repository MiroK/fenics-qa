from petsc4py import PETSc
from dolfin import *

mesh = UnitSquareMesh(100, 100)
V = FunctionSpace(mesh, 'CG', 2)
bc = DirichletBC(V, Constant(0.0), 'on_boundary')
dofmap = V.dofmap()
# Diagonal vector
v = interpolate(Expression('pow(x[0], 2)-x[1]', degree=2), V)
v = as_backend_type(v.vector()).vec()

comm = mesh.mpi_comm().tompi4py() 
mat = PETSc.Mat()
mat.create(comm)
# Local, global
sizes = [dofmap.index_map().size(IndexMap.MapSize_OWNED), 
         dofmap.index_map().size(IndexMap.MapSize_GLOBAL)]
# Square
mat.setSizes([sizes, sizes])
# Sparse
mat.setType('aij')
mat.setUp()
# Map from local rows to gloval rows
lgmap = map(int, dofmap.tabulate_local_to_global_dofs())
lgmap = PETSc.LGMap().create(lgmap, comm=comm)
mat.setLGMap(lgmap, lgmap)

# Fill the values
mat.setDiagonal(v)
mat.assemblyBegin()
mat.assemblyEnd()

A = PETScMatrix(mat)
bc.apply(A)

# Check if we can do matvec. Just no crash
x, y = mat.createVecs()
x.setRandom(); y.setRandom()
x, y = PETScVector(x), PETScVector(y)
A.mult(x, y)

# Verify the matrix
# The diagonal entry should be coordinates of the corresponing dof.
X = V.tabulate_dof_coordinates().reshape((-1, 2))
x, y = X[:, 0], X[:, 1]
values = x**2 - y  # Just like the interpolated expression
# For a row corresponding to a dof where bc was used the value should be 1
bc_dofs = bc.get_boundary_values().keys()
first, last = dofmap.ownership_range()
for ldof in range(last-first):
    gdof = dofmap.local_to_global_index(ldof)
    if first <= gdof < last:
        indices, row_values = A.getrow(gdof)
        assert len(indices) == len(row_values) == 1
        assert indices[0] == gdof

        if ldof in bc_dofs:
            assert int(row_values[0]) == 1
        else:
            assert near(values[ldof]-row_values[0], 0)
