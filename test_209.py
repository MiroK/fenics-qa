from dolfin import *

mesh = UnitSquareMesh(20, 20)

# Made up data for cell function
cell_f = CellFunction('size_t', mesh, 0)
for i, cell in enumerate(cells(mesh)): cell_f[cell] = i

V = FunctionSpace(mesh, 'DG', 0)
dofmap = V.dofmap()
# Represent data as function in V
f = Function(V)
vec = f.vector()

values = vec.get_local()
for cell in cells(mesh): 
    values[dofmap.cell_dofs(cell.index())[0]] = cell_f[cell]
vec.set_local(values)
vec.apply('insert')

plot(cell_f)
plot(f)
interactive()
