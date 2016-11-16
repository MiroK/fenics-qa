import dolfin as dlf

from petsc4py import PETSc

mesh = dlf.UnitSquareMesh(100, 100)
W = dlf.FunctionSpace(mesh, 'CG', 1)

def boundary(x, on_boundary):
    return on_boundary

bc = dlf.DirichletBC(W, dlf.Constant(0.0), boundary)

N = mesh.num_vertices()
v = PETSc.Vec()
v.create()
v.setSizes(N)
v.setType('standard')
v.setValues(range(N), [N]*N)

B_pet = PETSc.Mat()
B_pet.createAIJ([N,N], nnz=N)

lgmap = PETSc.LGMap().create(W.dofmap().dofs())
B_pet.setLGMap(lgmap, lgmap)

B_pet.setDiagonal(v)
B_pet.assemblyBegin()
B_pet.assemblyEnd()

B = dlf.PETScMatrix(B_pet)
bc.apply(B) # error
