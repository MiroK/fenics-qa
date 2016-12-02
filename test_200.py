from dolfin import *

mesh = Mesh()
editor = MeshEditor()
editor.open(mesh, 2, 2)

v = [0, 0, 0, 1, 0, 2, 1, 0, 1, 1, 1, 2, 2, 0, 2, 1, 2, 2]
num_vertices = len(v)/2
editor.init_vertices(num_vertices)

for i in range(num_vertices): editor.add_vertex(i, v[2*i], v[2*i+1])

num_cells = 8
editor.init_cells(num_cells);
editor.add_cell(0, 0,1,4);
editor.add_cell(1, 0,4,3);
editor.add_cell(2, 1,2,5);
editor.add_cell(3, 1,5,4);
editor.add_cell(4, 3,4,7);
editor.add_cell(5, 3,7,6);
editor.add_cell(6, 4,5,8);
editor.add_cell(7, 4,8,7);
editor.close();

V = VectorFunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = inner(Constant((1, 1)), v)*dx
bc = DirichletBC(V, Constant((0, 0)), 'on_boundary')
uh = Function(V)
solve(a == L, uh, bc)

